#!/usr/bin/env Rscript
# Convert RevPal.accdb -> revpal.sqlite and build a flat v_pasi table.
# Run from the project root:
#   Rscript app/convert.R

suppressPackageStartupMessages({
  library(DBI)
  library(odbc)
  library(RSQLite)
})

# Resolve the project root regardless of how the script is invoked.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_path <- if (length(file_arg)) file_arg else "app/convert.R"
here     <- normalizePath(dirname(script_path), mustWork = TRUE)
accdb    <- normalizePath(file.path(here, "..", "RevPal.accdb"), mustWork = TRUE)
sqlite_p <- file.path(here, "revpal.sqlite")

if (!file.exists(accdb)) stop("Cannot find ", accdb)
if (file.exists(sqlite_p)) file.remove(sqlite_p)

cat("Source:", accdb, "\n")
cat("Target:", sqlite_p, "\n\n")

src <- dbConnect(odbc::odbc(),
  .connection_string = sprintf(
    "Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;", accdb))
dst <- dbConnect(RSQLite::SQLite(), sqlite_p)

# --- 1. Copy real source tables only --------------------------------------
# The .accdb exposes ~150 objects: ~30 real `tbl*` tables of source data, and
# the rest are Access query views or scratch/temp tables. We whitelist `tbl*`
# and drop the obvious junk so the output stays small.
all_tabs <- dbListTables(src)
user_tabs <- all_tabs[
  grepl("^tbl", all_tabs) &
  !grepl("^tbl(ZTmp|Ztmp|Null)", all_tabs)
]

cat("Copying", length(user_tabs), "tables:\n")
for (t in user_tabs) {
  df <- tryCatch(dbReadTable(src, t), error = function(e) {
    cat(sprintf("  ! %-30s SKIPPED (%s)\n", t, conditionMessage(e))); NULL
  })
  if (is.null(df)) next
  dbWriteTable(dst, t, df, overwrite = TRUE)
  cat(sprintf("  - %-30s %6d rows\n", t, nrow(df)))
}

# --- 2. Build the flat v_pasi table ---------------------------------------
# PASI outcome IDs: 36 = PASI 50, 13 = PASI 75, 14 = PASI 90, 38 = PASI 100
# Study-char CatIDs: 49 = Cochrane Study ID, 57 = Drug, 58 = DoseAmount, 59 = DoseUnit

cat("\nBuilding v_pasi ...\n")

dbExecute(dst, "DROP TABLE IF EXISTS v_pasi")
dbExecute(dst, "
CREATE TABLE v_pasi AS
WITH
  trial AS (
    SELECT RefID, TextVal AS trial
    FROM   tblStudyChars
    WHERE  CatID = 49 AND ArmNo = 0
  ),
  drug AS (
    SELECT sc.RefID, sc.ArmNo, lk.AnsText AS drug
    FROM   tblStudyChars sc
    JOIN   tblStudyCharsDefsLookups lk
      ON   lk.CatID = 57 AND lk.AnsID = sc.ListVal
    WHERE  sc.CatID = 57 AND sc.ArmNo > 0
  ),
  dose_amount AS (
    SELECT RefID, ArmNo, NumVal AS dose_amount
    FROM   tblStudyChars
    WHERE  CatID = 58 AND ArmNo > 0
  ),
  dose_unit AS (
    SELECT sc.RefID, sc.ArmNo, lk.AnsText AS dose_unit
    FROM   tblStudyChars sc
    JOIN   tblStudyCharsDefsLookups lk
      ON   lk.CatID = 59 AND lk.AnsID = sc.ListVal
    WHERE  sc.CatID = 59 AND sc.ArmNo > 0
  ),
  pasi AS (
    -- For binary outcomes in this schema, responder count is in column k
    -- (not TP) and the arm total is in N.
    SELECT RefID, ArmNo, TimePeriod,
           MAX(CASE WHEN OutcomeID = 36 THEN k END) AS pasi50,
           MAX(CASE WHEN OutcomeID = 13 THEN k END) AS pasi75,
           MAX(CASE WHEN OutcomeID = 14 THEN k END) AS pasi90,
           MAX(CASE WHEN OutcomeID = 38 THEN k END) AS pasi100,
           MAX(N)                                   AS n
    FROM   tblIntraData
    WHERE  OutcomeID IN (36, 13, 14, 38)
      AND  COALESCE(SubgroupID, 0) = 0
    GROUP BY RefID, ArmNo, TimePeriod
  ),
  -- Time unit per RefID across the four PASI outcomes. In this dataset every
  -- ref uses one unit for all its PASI outcomes; if that ever stops being
  -- true MIN(Unit) just picks one deterministically. Default: weeks.
  tunit AS (
    SELECT ld.RefID, MIN(u.strUnit) AS unit
    FROM   tblLongitudinalDataDefs ld
    JOIN   tblLongitudinalUnitDefs  u ON u.UnitID = ld.Unit
    WHERE  ld.OutcomeID IN (36, 13, 14, 38)
    GROUP BY ld.RefID
  )
SELECT
  t.trial                                                       AS trial,
  p.RefID                                                       AS ref_id,
  p.ArmNo                                                       AS arm_no,
  a.ArmName                                                     AS arm_name,
  d.drug                                                        AS drug,
  -- Build dose like 45 mg or 0.5 mg/kg; strip trailing .0 from whole numbers.
  TRIM(
    CASE
      WHEN da.dose_amount IS NULL THEN ''
      WHEN da.dose_amount = CAST(da.dose_amount AS INTEGER)
        THEN CAST(CAST(da.dose_amount AS INTEGER) AS TEXT)
      ELSE CAST(da.dose_amount AS TEXT)
    END
    || ' ' || COALESCE(du.dose_unit, '')
  )                                                             AS dose,
  p.TimePeriod                                                  AS timepoint,
  COALESCE(tu.unit, 'wk')                                       AS timepoint_unit,
  p.n                                                           AS n,
  p.pasi50, p.pasi75, p.pasi90, p.pasi100
FROM       pasi p
LEFT JOIN  trial       t  ON t.RefID  = p.RefID
LEFT JOIN  tblArms     a  ON a.RefID  = p.RefID AND a.ArmNo  = p.ArmNo
LEFT JOIN  drug        d  ON d.RefID  = p.RefID AND d.ArmNo  = p.ArmNo
LEFT JOIN  dose_amount da ON da.RefID = p.RefID AND da.ArmNo = p.ArmNo
LEFT JOIN  dose_unit   du ON du.RefID = p.RefID AND du.ArmNo = p.ArmNo
LEFT JOIN  tunit       tu ON tu.RefID = p.RefID
ORDER BY   trial, p.ArmNo, p.TimePeriod
")

n_rows <- dbGetQuery(dst, "SELECT COUNT(*) AS n FROM v_pasi")$n
n_drug <- dbGetQuery(dst, "SELECT COUNT(DISTINCT drug) AS n FROM v_pasi WHERE drug IS NOT NULL")$n
cat(sprintf("  v_pasi rows: %d\n  distinct drugs: %d\n", n_rows, n_drug))

cat("\nSample (first 8 rows):\n")
print(dbGetQuery(dst, "SELECT trial, drug, dose, timepoint, n, pasi50, pasi75, pasi90, pasi100
                       FROM v_pasi LIMIT 8"))

dbDisconnect(src)

# Reclaim space and compact the file.
dbExecute(dst, "VACUUM")
dbDisconnect(dst)

cat(sprintf("\nDone. SQLite file: %s (%.1f MB)\n",
            sqlite_p, file.info(sqlite_p)$size / 1024 / 1024))
