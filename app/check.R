#!/usr/bin/env Rscript
# Data-quality checks for revpal.sqlite. Run from the project root:
#   Rscript app/check.R
#
# Prints each problem class with a short heading and the offending rows.
# Empty sections are summarised as "OK".

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_path <- if (length(file_arg)) file_arg else "app/check.R"
here     <- normalizePath(dirname(script_path), mustWork = TRUE)
sqlite_p <- file.path(here, "revpal.sqlite")
if (!file.exists(sqlite_p)) stop("revpal.sqlite not found - run convert.R first.")

con <- dbConnect(SQLite(), sqlite_p, flags = SQLITE_RO)
on.exit(dbDisconnect(con), add = TRUE)

options(width = 200)

q <- function(sql) dbGetQuery(con, sql)

# Truncate long character columns so wide arm_name fields don't blow up the
# console layout. Numeric/integer columns are left alone.
trunc_chars <- function(df, max_chars = 60) {
  for (nm in names(df)) {
    if (is.character(df[[nm]])) {
      long <- !is.na(df[[nm]]) & nchar(df[[nm]]) > max_chars
      df[[nm]][long] <- paste0(substr(df[[nm]][long], 1, max_chars - 1), "...")
    }
  }
  df
}

# Pretty-print one check. `df` is the offending rows (0 rows => OK).
# `max_rows` caps the printed rows for noisy informational checks.
section <- function(title, df, note = NULL, max_rows = Inf) {
  cat("\n--- ", title, " ---\n", sep = "")
  if (!is.null(note)) cat(note, "\n", sep = "")
  if (nrow(df) == 0) {
    cat("OK (0 problems)\n")
  } else {
    cat(sprintf("%d problem row(s):\n", nrow(df)))
    shown <- trunc_chars(if (nrow(df) > max_rows) head(df, max_rows) else df)
    print(shown, row.names = FALSE, right = FALSE)
    if (nrow(df) > max_rows) {
      cat(sprintf("  ... %d more rows not shown.\n", nrow(df) - max_rows))
    }
  }
}

cat("Checking", sqlite_p, "\n")

# 1. Arm has dose or PASI data but no drug name.
section(
  "Missing drug (arm has dose or PASI data)",
  q("
    SELECT DISTINCT trial, ref_id, arm_no, arm_name, dose
    FROM   v_pasi
    WHERE  (drug IS NULL OR drug = '')
      AND  (
        (dose IS NOT NULL AND dose != '')
        OR pasi50  IS NOT NULL OR pasi75  IS NOT NULL
        OR pasi90  IS NOT NULL OR pasi100 IS NOT NULL
      )
    ORDER BY trial, ref_id, arm_no
  ")
)

# 2. PASI monotonicity. More-lenient thresholds must have >= responders than
#    stricter ones: pasi50 >= pasi75 >= pasi90 >= pasi100.
section(
  "PASI monotonicity violations (e.g. PASI 90 > PASI 75)",
  q("
    SELECT trial, ref_id, arm_no, drug, dose, timepoint,
           pasi50, pasi75, pasi90, pasi100
    FROM   v_pasi
    WHERE  (pasi75  IS NOT NULL AND pasi50  IS NOT NULL AND pasi75  > pasi50)
       OR  (pasi90  IS NOT NULL AND pasi75  IS NOT NULL AND pasi90  > pasi75)
       OR  (pasi100 IS NOT NULL AND pasi90  IS NOT NULL AND pasi100 > pasi90)
       OR  (pasi90  IS NOT NULL AND pasi50  IS NOT NULL AND pasi90  > pasi50)
       OR  (pasi100 IS NOT NULL AND pasi75  IS NOT NULL AND pasi100 > pasi75)
       OR  (pasi100 IS NOT NULL AND pasi50  IS NOT NULL AND pasi100 > pasi50)
    ORDER BY trial, ref_id, arm_no, timepoint
  ")
)

# 3. Responder count exceeds arm N.
section(
  "Responder count exceeds N",
  q("
    SELECT trial, ref_id, arm_no, drug, timepoint, n,
           pasi50, pasi75, pasi90, pasi100
    FROM   v_pasi
    WHERE  n IS NOT NULL AND n > 0
      AND  (pasi50  > n OR pasi75 > n OR pasi90 > n OR pasi100 > n)
    ORDER BY trial, ref_id, arm_no, timepoint
  ")
)

# 4. PASI responder count present but N missing or zero (can't compute %).
section(
  "PASI responders present but N missing/zero",
  q("
    SELECT trial, ref_id, arm_no, drug, timepoint, n,
           pasi50, pasi75, pasi90, pasi100
    FROM   v_pasi
    WHERE  (n IS NULL OR n = 0)
      AND  (pasi50 IS NOT NULL OR pasi75 IS NOT NULL
            OR pasi90 IS NOT NULL OR pasi100 IS NOT NULL)
    ORDER BY trial, ref_id, arm_no, timepoint
  ")
)

# 5. Inconsistent N across the four PASI outcomes within one arm × timepoint.
#    v_pasi takes MAX(N); if the per-outcome N values disagree, the curator
#    likely typo'd one of them.
section(
  "N differs between PASI outcomes at the same arm × timepoint",
  q("
    SELECT RefID AS ref_id, ArmNo AS arm_no, TimePeriod AS timepoint,
           MIN(N) AS n_min, MAX(N) AS n_max, COUNT(DISTINCT N) AS distinct_n
    FROM   tblIntraData
    WHERE  OutcomeID IN (36, 13, 14, 38) AND N IS NOT NULL
    GROUP  BY RefID, ArmNo, TimePeriod
    HAVING COUNT(DISTINCT N) > 1
    ORDER BY ref_id, arm_no, timepoint
  "),
  note = "(v_pasi reports MAX(N); these arms have conflicting Ns in the raw data.)"
)

# 6. Missing trial (Cochrane Study ID) for an arm that has data.
section(
  "Missing Cochrane Study ID",
  q("
    SELECT DISTINCT ref_id, arm_no, arm_name, drug
    FROM   v_pasi
    WHERE  trial IS NULL OR trial = ''
    ORDER BY ref_id, arm_no
  ")
)

# 7. Missing timepoint.
section(
  "Missing timepoint",
  q("
    SELECT trial, ref_id, arm_no, drug, timepoint, n,
           pasi50, pasi75, pasi90, pasi100
    FROM   v_pasi
    WHERE  timepoint IS NULL
    ORDER BY trial, ref_id, arm_no
  ")
)

# 8. Dose amount with no unit, or unit with no amount.
section(
  "Dose amount/unit mismatch",
  q("
    WITH
      da AS (SELECT RefID, ArmNo, NumVal AS amt
             FROM tblStudyChars WHERE CatID = 58 AND ArmNo > 0),
      du AS (SELECT RefID, ArmNo, ListVal AS unit
             FROM tblStudyChars WHERE CatID = 59 AND ArmNo > 0)
    SELECT COALESCE(da.RefID, du.RefID) AS ref_id,
           COALESCE(da.ArmNo, du.ArmNo) AS arm_no,
           da.amt  AS dose_amount,
           du.unit AS dose_unit_id
    FROM   da LEFT JOIN du ON da.RefID = du.RefID AND da.ArmNo = du.ArmNo
    WHERE  da.amt IS NOT NULL AND du.unit IS NULL
    UNION ALL
    SELECT COALESCE(da.RefID, du.RefID) AS ref_id,
           COALESCE(da.ArmNo, du.ArmNo) AS arm_no,
           da.amt  AS dose_amount,
           du.unit AS dose_unit_id
    FROM   du LEFT JOIN da ON da.RefID = du.RefID AND da.ArmNo = du.ArmNo
    WHERE  du.unit IS NOT NULL AND da.amt IS NULL
    ORDER BY ref_id, arm_no
  ")
)

# 9. Duplicate (ref_id, arm_no, timepoint) in v_pasi — should be impossible
#    given the GROUP BY, but worth confirming.
section(
  "Duplicate (ref_id, arm_no, timepoint) rows in v_pasi",
  q("
    SELECT ref_id, arm_no, timepoint, COUNT(*) AS n_rows
    FROM   v_pasi
    GROUP  BY ref_id, arm_no, timepoint
    HAVING COUNT(*) > 1
    ORDER BY ref_id, arm_no, timepoint
  ")
)

# 10. Same drug name spelled differently (whitespace / case). Informational —
#     might be legitimate (e.g. brand vs INN) but often a typo.
section(
  "Drug names that collapse under trim+lower (possible typos)",
  q("
    WITH norm AS (
      SELECT drug,
             LOWER(TRIM(drug)) AS key
      FROM   v_pasi
      WHERE  drug IS NOT NULL AND drug != ''
      GROUP  BY drug
    )
    SELECT key, GROUP_CONCAT(DISTINCT drug) AS variants, COUNT(*) AS n_variants
    FROM   norm
    GROUP  BY key
    HAVING COUNT(*) > 1
    ORDER BY key
  ")
)

# 11. Arms in tblArms with no PASI data at all. Informational, not an error
#     (some trials extract other outcomes), but useful to scan.
section(
  "Arms with no PASI data (informational)",
  q("
    SELECT a.RefID AS ref_id, a.ArmNo AS arm_no, a.ArmName AS arm_name
    FROM   tblArms a
    LEFT   JOIN v_pasi v ON v.ref_id = a.RefID AND v.arm_no = a.ArmNo
    WHERE  v.ref_id IS NULL
    ORDER BY ref_id, arm_no
  "),
  max_rows = 20
)

cat("\nDone.\n")
