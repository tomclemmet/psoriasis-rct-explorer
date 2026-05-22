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

# --- 2. Build the flat per-endpoint-group views ---------------------------
# Each view is one row per (RefID, ArmNo, TimePeriod), pivoted wide on the
# outcomes in that group. Common context (trial / drug / dose / timepoint
# unit) is shared via the CTEs below; each view supplies its own `outcomes`
# CTE that picks and pivots the relevant OutcomeIDs.
#
# Study-char CatIDs: 49 = Cochrane Study ID, 57 = Drug, 58 = DoseAmount,
#                    59 = DoseUnit.
# Outcome IDs are listed alongside each view's `outcomes` CTE.

context_ctes <- "
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
  )
"

# Time-unit CTE — parameterised by which OutcomeIDs to consult. Same unit
# across outcomes within a trial in this data; MIN(Unit) is deterministic.
tunit_cte <- function(outcome_ids) {
  ids <- paste(outcome_ids, collapse = ", ")
  sprintf("
  tunit AS (
    SELECT ld.RefID, MIN(u.strUnit) AS unit
    FROM   tblLongitudinalDataDefs ld
    JOIN   tblLongitudinalUnitDefs  u ON u.UnitID = ld.Unit
    WHERE  ld.OutcomeID IN (%s)
    GROUP BY ld.RefID
  )", ids)
}

# Standard SELECT prefix (trial, ref_id, arm_no, arm_name, drug, dose,
# timepoint, timepoint_unit, n). The caller appends the view-specific columns.
common_select_head <- "
SELECT
  t.trial                                                       AS trial,
  o.RefID                                                       AS ref_id,
  o.ArmNo                                                       AS arm_no,
  a.ArmName                                                     AS arm_name,
  d.drug                                                        AS drug,
  TRIM(
    CASE
      WHEN da.dose_amount IS NULL THEN ''
      WHEN da.dose_amount = CAST(da.dose_amount AS INTEGER)
        THEN CAST(CAST(da.dose_amount AS INTEGER) AS TEXT)
      ELSE CAST(da.dose_amount AS TEXT)
    END
    || ' ' || COALESCE(du.dose_unit, '')
  )                                                             AS dose,
  o.TimePeriod                                                  AS timepoint,
  COALESCE(tu.unit, 'wk')                                       AS timepoint_unit,
  o.n                                                           AS n
"

common_join_tail <- "
FROM       outcomes o
LEFT JOIN  trial       t  ON t.RefID  = o.RefID
LEFT JOIN  tblArms     a  ON a.RefID  = o.RefID AND a.ArmNo  = o.ArmNo
LEFT JOIN  drug        d  ON d.RefID  = o.RefID AND d.ArmNo  = o.ArmNo
LEFT JOIN  dose_amount da ON da.RefID = o.RefID AND da.ArmNo = o.ArmNo
LEFT JOIN  dose_unit   du ON du.RefID = o.RefID AND du.ArmNo = o.ArmNo
LEFT JOIN  tunit       tu ON tu.RefID = o.RefID
ORDER BY   trial, o.ArmNo, o.TimePeriod
"

# Build and run one view.
#   name        - destination table name in SQLite
#   outcome_ids - OutcomeIDs that contribute (used for the tunit CTE)
#   outcomes_sql - the `outcomes` CTE body (everything between
#                  `outcomes AS (` and the closing `)`). Must produce columns
#                  RefID, ArmNo, TimePeriod, n, and whatever view-specific
#                  pivot columns the SELECT references.
#   extra_cols  - the comma-leading list of `o.colname` expressions appended
#                 to common_select_head.
build_view <- function(name, outcome_ids, outcomes_sql, extra_cols) {
  cat(sprintf("\nBuilding %s ...\n", name))
  dbExecute(dst, sprintf("DROP TABLE IF EXISTS %s", name))
  sql <- paste0(
    sprintf("CREATE TABLE %s AS\nWITH", name),
    context_ctes, ",",
    tunit_cte(outcome_ids), ",\n",
    "  outcomes AS (\n", outcomes_sql, "\n  )\n",
    common_select_head,
    ",", extra_cols, "\n",
    common_join_tail
  )
  dbExecute(dst, sql)
  n_rows <- dbGetQuery(dst, sprintf("SELECT COUNT(*) AS n FROM %s", name))$n
  cat(sprintf("  %s rows: %d\n", name, n_rows))
}

# v_pasi — binary PASI thresholds (existing).
build_view(
  name = "v_pasi",
  outcome_ids = c(36, 13, 14, 38),
  outcomes_sql = "
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
  ",
  extra_cols = "o.pasi50, o.pasi75, o.pasi90, o.pasi100"
)

# v_pasi_abs — continuous PASI (Absolute PASI, change from baseline).
# 46 = Absolute PASI, 34 = Absolute PASI change from baseline.
build_view(
  name = "v_pasi_abs",
  outcome_ids = c(46, 34),
  outcomes_sql = "
    SELECT RefID, ArmNo, TimePeriod,
           MAX(CASE WHEN OutcomeID = 46 THEN Mean    END) AS abs_pasi_mean,
           MAX(CASE WHEN OutcomeID = 46 THEN SD      END) AS abs_pasi_sd,
           MAX(CASE WHEN OutcomeID = 46 THEN Median  END) AS abs_pasi_median,
           MAX(CASE WHEN OutcomeID = 46 THEN loIQR   END) AS abs_pasi_lo_iqr,
           MAX(CASE WHEN OutcomeID = 46 THEN hiIQR   END) AS abs_pasi_hi_iqr,
           MAX(CASE WHEN OutcomeID = 34 THEN Mean    END) AS abs_pasi_change_mean,
           MAX(CASE WHEN OutcomeID = 34 THEN SD      END) AS abs_pasi_change_sd,
           MAX(CASE WHEN OutcomeID = 34 THEN Median  END) AS abs_pasi_change_median,
           MAX(N)                                          AS n
    FROM   tblIntraData
    WHERE  OutcomeID IN (46, 34)
      AND  COALESCE(SubgroupID, 0) = 0
    GROUP BY RefID, ArmNo, TimePeriod
  ",
  extra_cols = "
  o.abs_pasi_mean, o.abs_pasi_sd, o.abs_pasi_median,
  o.abs_pasi_lo_iqr, o.abs_pasi_hi_iqr,
  o.abs_pasi_change_mean, o.abs_pasi_change_sd, o.abs_pasi_change_median"
)

# v_dlqi — mixed binary + continuous DLQI endpoints.
# Binary: 41 = DLQI 0 or 1, 51 = DLQI 0, 35 = 5+ pt decrease,
#         50 = 4+ pt decrease, 112 = DLQI <= 5 (5+ at baseline).
# Continuous: 43 = Absolute DLQI, 56 = Absolute change in DLQI.
build_view(
  name = "v_dlqi",
  outcome_ids = c(41, 51, 35, 50, 112, 43, 56),
  outcomes_sql = "
    SELECT RefID, ArmNo, TimePeriod,
           MAX(CASE WHEN OutcomeID = 41  THEN k END) AS dlqi_0_1,
           MAX(CASE WHEN OutcomeID = 51  THEN k END) AS dlqi_0,
           MAX(CASE WHEN OutcomeID = 35  THEN k END) AS dlqi_5pt_dec,
           MAX(CASE WHEN OutcomeID = 50  THEN k END) AS dlqi_4pt_dec,
           MAX(CASE WHEN OutcomeID = 112 THEN k END) AS dlqi_le5,
           MAX(CASE WHEN OutcomeID = 43  THEN Mean   END) AS abs_dlqi_mean,
           MAX(CASE WHEN OutcomeID = 43  THEN SD     END) AS abs_dlqi_sd,
           MAX(CASE WHEN OutcomeID = 43  THEN Median END) AS abs_dlqi_median,
           MAX(CASE WHEN OutcomeID = 56  THEN Mean   END) AS abs_dlqi_change_mean,
           MAX(CASE WHEN OutcomeID = 56  THEN SD     END) AS abs_dlqi_change_sd,
           MAX(CASE WHEN OutcomeID = 56  THEN Median END) AS abs_dlqi_change_median,
           MAX(N)                                        AS n
    FROM   tblIntraData
    WHERE  OutcomeID IN (41, 51, 35, 50, 112, 43, 56)
      AND  COALESCE(SubgroupID, 0) = 0
    GROUP BY RefID, ArmNo, TimePeriod
  ",
  extra_cols = "
  o.dlqi_0_1, o.dlqi_0, o.dlqi_5pt_dec, o.dlqi_4pt_dec, o.dlqi_le5,
  o.abs_dlqi_mean, o.abs_dlqi_sd, o.abs_dlqi_median,
  o.abs_dlqi_change_mean, o.abs_dlqi_change_sd, o.abs_dlqi_change_median"
)

# v_safety — common binary safety outcomes.
# 20 = SAE, 48 = Disc (any), 37 = Disc (AE), 23 = Serious infection,
# 21 = Injection site reaction, 96 = Malignancy, 24 = NMSC,
# 25 = Malignancy non-NMSC.
build_view(
  name = "v_safety",
  outcome_ids = c(20, 48, 37, 23, 21, 96, 24, 25),
  outcomes_sql = "
    SELECT RefID, ArmNo, TimePeriod,
           MAX(CASE WHEN OutcomeID = 20 THEN k END) AS sae,
           MAX(CASE WHEN OutcomeID = 48 THEN k END) AS disc_any,
           MAX(CASE WHEN OutcomeID = 37 THEN k END) AS disc_ae,
           MAX(CASE WHEN OutcomeID = 23 THEN k END) AS serious_infection,
           MAX(CASE WHEN OutcomeID = 21 THEN k END) AS injection_site_rxn,
           MAX(CASE WHEN OutcomeID = 96 THEN k END) AS malignancy,
           MAX(CASE WHEN OutcomeID = 24 THEN k END) AS nmsc,
           MAX(CASE WHEN OutcomeID = 25 THEN k END) AS malignancy_non_nmsc,
           MAX(N)                                    AS n
    FROM   tblIntraData
    WHERE  OutcomeID IN (20, 48, 37, 23, 21, 96, 24, 25)
      AND  COALESCE(SubgroupID, 0) = 0
    GROUP BY RefID, ArmNo, TimePeriod
  ",
  extra_cols = "
  o.sae, o.disc_any, o.disc_ae, o.serious_infection,
  o.injection_site_rxn, o.malignancy, o.nmsc, o.malignancy_non_nmsc"
)

n_drug <- dbGetQuery(dst, "SELECT COUNT(DISTINCT drug) AS n FROM v_pasi WHERE drug IS NOT NULL")$n
cat(sprintf("\nDistinct drugs in v_pasi: %d\n", n_drug))

cat("\nSample v_pasi (first 5 rows):\n")
print(dbGetQuery(dst, "SELECT trial, drug, dose, timepoint, n, pasi50, pasi75, pasi90, pasi100
                       FROM v_pasi LIMIT 5"))

dbDisconnect(src)

# Reclaim space and compact the file.
dbExecute(dst, "VACUUM")
dbDisconnect(dst)

cat(sprintf("\nDone. SQLite file: %s (%.1f MB)\n",
            sqlite_p, file.info(sqlite_p)$size / 1024 / 1024))
