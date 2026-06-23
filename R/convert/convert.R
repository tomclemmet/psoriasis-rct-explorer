#!/usr/bin/env Rscript
# Convert RevPal.accdb -> psoriasis-rcts.sqlite as a normalised relational
# database, then expose the wide SQL views the Shiny app reads (built in
# views.R). Run from RStudio (open in the project) or the terminal:
#   Rscript app/convert.R
#
# Output schema (one fact per row, no repeated text):
#
#   Lookups
#     drugs(drug_id PK, drug_name UNIQUE)
#     dose_units(unit_id PK, unit_name UNIQUE)
#     data_types(data_type_id PK, name)            -- from tblDataTypes
#     outcomes(outcome_id PK, code, label, subcategory,
#              data_type_id FK, endpoint_group)    -- full Access catalogue
#     subgroups(subgroup_id PK, subgroup_name)     -- only the "Main analysis"
#                                                  -- sentinel (0) is kept
#
#   Entities
#     studies(study_id PK, trial, timepoint_unit)
#         -- Keyed by the *primary* publication's RefID. timepoint_unit is the
#         -- per-study unit ("wk" almost always) for measurement timepoints.
#     publications(publication_id PK, study_id FK, is_primary, doi, title,
#                  authors, year, journal, volume, issue, page_start, page_end,
#                  notes)
#     arms(arm_id PK, study_id FK, arm_no, arm_name,
#          drug_id FK, dose_amount, dose_unit_id FK)
#
#   Facts
#     measurements(measurement_id PK, arm_id FK, outcome_id FK,
#                  subgroup_id FK, timepoint, k, n, mean, sd, median,
#                  lo_iqr, hi_iqr)
#         -- The 21 wide-view outcomes plus baseline PASI (outcome 11), main
#         -- analysis only, timepoints within 24 week-equivalent weeks.
#
#   Views (views.R): v_pasi, v_dlqi, v_safety. Each carries its tab's binary
#   responder columns plus any absolute / change-from-baseline columns.
#
# The Access file is the source of truth and is never modified. This script is
# slimmed to exactly what the app reads now: no subgroup rows/tables, no
# demographic/baseline-characteristic rows (except baseline PASI), no
# descriptive trial-info columns, and timepoints capped at 24 weeks.

suppressPackageStartupMessages({
  library(DBI)
  library(odbc)
  library(RSQLite)
})

# Drop measurements beyond this many week-equivalent weeks.
MAX_TIMEPOINT_WK <- 24
# Baseline PASI is a "Psoriasis characteristics" outcome; always kept (the app
# uses it as the Absolute-PASI baseline) even though it has no view `code`.
BASELINE_PASI_OUTCOME_ID <- 11L

accdb    <- normalizePath("RevPal.accdb", mustWork = TRUE)
sqlite_p <- "app/psoriasis-rcts.sqlite"

if (file.exists(sqlite_p) && !file.remove(sqlite_p)) {
  stop("Could not remove existing ", sqlite_p,
       " - is the Shiny app or another R session holding it open?")
}

cat("Source:", accdb, "\n")
cat("Target:", sqlite_p, "\n\n")

src <- dbConnect(odbc::odbc(),
  .connection_string = sprintf(
    "Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;", accdb))
dst <- dbConnect(RSQLite::SQLite(), sqlite_p)

dbExecute(dst, "PRAGMA foreign_keys = ON")

# --- 1. Outcome catalogue (wide-view decoration) ---------------------------
# Outcome IDs match the Access OutcomeID. `code` is the short identifier used
# inside the views (and exposed as the view's column name); `endpoint_group`
# names the tab. Only the 21 outcomes used by a view appear here; the full
# catalogue is read from tblOutcomeDefs below and merged with this decoration.
wide_outcomes_df <- read.table(text = "
outcome_id|code|endpoint_group
36|pasi50|pasi
13|pasi75|pasi
14|pasi90|pasi
38|pasi100|pasi
46|abs_pasi|pasi
34|abs_pasi_change|pasi
41|dlqi_0_1|dlqi
51|dlqi_0|dlqi
35|dlqi_5pt_dec|dlqi
50|dlqi_4pt_dec|dlqi
112|dlqi_le5|dlqi
43|abs_dlqi|dlqi
56|abs_dlqi_change|dlqi
20|sae|safety
48|disc_any|safety
37|disc_ae|safety
23|serious_infection|safety
21|injection_site_rxn|safety
96|malignancy|safety
24|nmsc|safety
25|malignancy_non_nmsc|safety
", sep = "|", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)

# --- 2. Read what we need from Access --------------------------------------
cat("Reading source tables...\n")
intra   <- dbReadTable(src, "tblIntraData")
chars   <- dbReadTable(src, "tblStudyChars")
lookups <- dbReadTable(src, "tblStudyCharsDefsLookups")
arms_x  <- dbReadTable(src, "tblArms")
refs    <- dbReadTable(src, "tblRefs")
lddefs  <- dbReadTable(src, "tblLongitudinalDataDefs")
lunits  <- dbReadTable(src, "tblLongitudinalUnitDefs")
odefs   <- dbReadTable(src, "tblOutcomeDefs")    # full outcomes catalogue
dtypes  <- dbReadTable(src, "tblDataTypes")      # DataTypeID -> name
studefs <- dbReadTable(src, "tblStudyDefs")      # RefID, ParentID, IncExc, ...
dbDisconnect(src)

# --- 2b. Curator exclusion flags ------------------------------------------
# tblRefs.Duplicate flags publication rows that are duplicate citations of
# another row for the same trial. Drop these so they never appear.
refs <- refs[is.na(refs$Duplicate) | refs$Duplicate != 1, ]

# tblStudyDefs.IncExc:
#   NULL / 0 -> include normally
#   1        -> excluded as secondary publication (captured via ParentID in
#               `publications`, not in `studies`)
#   2 - 4    -> excluded for other reasons. Drop the study and every dependent
#               row.
excluded_study_ids <- studefs$RefID[
  (is.na(studefs$ParentID) | studefs$ParentID == 0) &
  !is.na(studefs$IncExc) & studefs$IncExc %in% c(2L, 3L, 4L)
]

# --- 3. Create destination schema -----------------------------------------
ddl <- "
CREATE TABLE drugs (
  drug_id    INTEGER PRIMARY KEY,
  drug_name  TEXT NOT NULL UNIQUE
);
CREATE TABLE dose_units (
  unit_id    INTEGER PRIMARY KEY,
  unit_name  TEXT NOT NULL UNIQUE
);
CREATE TABLE data_types (
  data_type_id  INTEGER PRIMARY KEY,
  name          TEXT NOT NULL UNIQUE
);
CREATE TABLE outcomes (
  outcome_id      INTEGER PRIMARY KEY,
  code            TEXT UNIQUE,           -- short identifier for the 21 wide-view
                                         -- outcomes; NULL for baseline chars etc.
  label           TEXT NOT NULL,         -- OutcomeName from tblOutcomeDefs
  subcategory     TEXT,                  -- 'Demographics', 'PASI', 'DLQI', ...
  data_type_id    INTEGER REFERENCES data_types(data_type_id),
  endpoint_group  TEXT                   -- 'pasi','dlqi','safety', or NULL
);
CREATE TABLE subgroups (
  subgroup_id    INTEGER PRIMARY KEY,
  subgroup_name  TEXT
);
CREATE TABLE studies (
  study_id        INTEGER PRIMARY KEY,   -- = primary publication's RefID
  trial           TEXT,                  -- Cochrane Study ID (CatID 49)
  timepoint_unit  TEXT                   -- per-study timepoint unit ('wk' etc.)
);
CREATE TABLE publications (
  publication_id  INTEGER PRIMARY KEY,         -- = tblRefs.ID
  study_id        INTEGER NOT NULL REFERENCES studies(study_id),
  is_primary      INTEGER NOT NULL,            -- 1 for the primary pub, else 0
  doi             TEXT,
  title           TEXT,
  authors         TEXT,
  year            INTEGER,                     -- parsed from tblRefs.Date
  journal         TEXT,
  volume          TEXT,                        -- tblRefs.Vol (occasionally non-numeric)
  issue           TEXT,
  page_start      TEXT,                        -- e.g. S12
  page_end        TEXT,
  notes           TEXT
);
CREATE INDEX idx_publications_study ON publications(study_id);
CREATE TABLE arms (
  arm_id        INTEGER PRIMARY KEY,
  study_id      INTEGER NOT NULL REFERENCES studies(study_id),
  arm_no        INTEGER NOT NULL,
  arm_name      TEXT,
  drug_id       INTEGER REFERENCES drugs(drug_id),
  dose_amount   REAL,
  dose_unit_id  INTEGER REFERENCES dose_units(unit_id),
  UNIQUE (study_id, arm_no)
);
CREATE TABLE measurements (
  measurement_id  INTEGER PRIMARY KEY,
  arm_id          INTEGER NOT NULL REFERENCES arms(arm_id),
  outcome_id      INTEGER NOT NULL REFERENCES outcomes(outcome_id),
  subgroup_id     INTEGER NOT NULL REFERENCES subgroups(subgroup_id),
  timepoint       INTEGER,
  k               INTEGER,
  n               INTEGER,
  mean            REAL,
  sd              REAL,
  median          REAL,
  lo_iqr          REAL,
  hi_iqr          REAL,
  UNIQUE (arm_id, outcome_id, subgroup_id, timepoint)
);
CREATE INDEX idx_measurements_outcome_arm ON measurements(outcome_id, arm_id);
"
for (stmt in strsplit(ddl, ";\\s*\n", perl = TRUE)[[1]]) {
  s <- trimws(stmt)
  if (nzchar(s)) dbExecute(dst, s)
}

# --- 4. Seed data_types, outcomes, subgroups ------------------------------
data_types_df <- data.frame(
  data_type_id = dtypes$DataTypeID,
  name         = dtypes$DataType,
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "data_types", data_types_df, append = TRUE)

# Full outcome catalogue from tblOutcomeDefs, decorated with view metadata for
# the 21 outcomes the app pivots into v_pasi / v_dlqi / v_safety.
outcomes_df <- merge(
  data.frame(outcome_id   = odefs$OutcomeID,
             label        = odefs$OutcomeName,
             subcategory  = odefs$SubCategory,
             data_type_id = odefs$DataTypeID,
             stringsAsFactors = FALSE),
  wide_outcomes_df,
  by = "outcome_id", all.x = TRUE
)
outcomes_df <- outcomes_df[, c("outcome_id","code","label","subcategory",
                               "data_type_id","endpoint_group")]
outcomes_df <- outcomes_df[order(outcomes_df$outcome_id), ]
dbWriteTable(dst, "outcomes", outcomes_df, append = TRUE)

# Subgroups: only the main-analysis sentinel (0). Every measurement carries
# subgroup_id = 0.
dbWriteTable(dst, "subgroups",
             data.frame(subgroup_id = 0L, subgroup_name = "Main analysis",
                        stringsAsFactors = FALSE),
             append = TRUE)

# --- 5. Lookup population -------------------------------------------------
# Drugs: distinct lookup labels (CatID = 57) actually referenced by an arm.
drug_used  <- chars[chars$CatID == 57 & !is.na(chars$ArmNo) & chars$ArmNo > 0, ]
drug_names <- sort(unique(lookups$AnsText[
  lookups$CatID == 57 & lookups$AnsID %in% drug_used$ListVal
]))
drug_names <- drug_names[!is.na(drug_names) & nzchar(drug_names)]
drugs_df <- data.frame(
  drug_id   = seq_along(drug_names),
  drug_name = drug_names,
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "drugs", drugs_df, append = TRUE)
drug_id_of <- setNames(drugs_df$drug_id, drugs_df$drug_name)

# Dose units (CatID = 59).
unit_used  <- chars[chars$CatID == 59 & !is.na(chars$ArmNo) & chars$ArmNo > 0, ]
unit_names <- sort(unique(lookups$AnsText[
  lookups$CatID == 59 & lookups$AnsID %in% unit_used$ListVal
]))
unit_names <- unit_names[!is.na(unit_names) & nzchar(unit_names)]
dose_units_df <- data.frame(
  unit_id   = seq_along(unit_names),
  unit_name = unit_names,
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "dose_units", dose_units_df, append = TRUE)
dose_unit_id_of <- setNames(dose_units_df$unit_id, dose_units_df$unit_name)

# --- 6. Studies -----------------------------------------------------------
# Universe = *primary* RefIDs (ParentID NULL/0) that have any extracted data,
# minus the IncExc-excluded ones. Secondaries never get a `studies` row; their
# bibliography lives in `publications`.
all_primary_ids <- studefs$RefID[is.na(studefs$ParentID) | studefs$ParentID == 0]
primary_ids <- setdiff(all_primary_ids, excluded_study_ids)

data_ref_ids <- unique(as.integer(c(
  intra$RefID,
  chars$RefID[chars$RefID %in% primary_ids],
  arms_x$RefID
)))
data_ref_ids <- data_ref_ids[!is.na(data_ref_ids)]
all_ref_ids <- sort(intersect(primary_ids, data_ref_ids))

# Trial name (CatID 49): study-level TextVal (ArmNo = 0), keyed by RefID.
trial_rows <- chars[chars$CatID == 49 & (is.na(chars$ArmNo) | chars$ArmNo == 0), ]
trial_of   <- setNames(trial_rows$TextVal, as.character(trial_rows$RefID))

# Per-study timepoint unit: MIN(strUnit) across every outcome present in
# tblLongitudinalDataDefs for the study (deterministic, constant per study
# here). Defaults to "wk". Used both for the studies.timepoint_unit column and
# the week-equivalent cutoff applied to measurements below.
ld_join <- merge(
  lddefs[, c("RefID", "Unit")],
  lunits[, c("UnitID", "strUnit")],
  by.x = "Unit", by.y = "UnitID", all.x = TRUE
)
tp_unit_by_ref <- tapply(ld_join$strUnit, ld_join$RefID,
                         function(x) min(x[!is.na(x)]))
study_tp_unit_name <- setNames(
  ifelse(is.na(tp_unit_by_ref[as.character(all_ref_ids)]),
         "wk", tp_unit_by_ref[as.character(all_ref_ids)]),
  as.character(all_ref_ids)
)

studies_df <- data.frame(
  study_id       = all_ref_ids,
  trial          = unname(trial_of[as.character(all_ref_ids)]),
  timepoint_unit = unname(study_tp_unit_name[as.character(all_ref_ids)]),
  stringsAsFactors = FALSE
)
# Blank trial strings -> NA.
studies_df$trial[!is.na(studies_df$trial) & !nzchar(studies_df$trial)] <- NA
dbWriteTable(dst, "studies", studies_df, append = TRUE)

# --- 6b. Publications -----------------------------------------------------
# One row per tblRefs entry, pointing at its primary's study_id (= ParentID if
# non-NULL/>0, else its own RefID). Only keep publications whose resolved
# primary has a `studies` row.
parent_of <- setNames(
  ifelse(is.na(studefs$ParentID) | studefs$ParentID == 0,
         studefs$RefID, studefs$ParentID),
  as.character(studefs$RefID)
)
pub_study_id <- unname(parent_of[as.character(refs$ID)])
pub_study_id[is.na(pub_study_id)] <- refs$ID[is.na(pub_study_id)]

publications_df <- data.frame(
  publication_id = as.integer(refs$ID),
  study_id       = as.integer(pub_study_id),
  is_primary     = as.integer(refs$ID == pub_study_id),
  doi            = refs$DOI,
  title          = refs$Title,
  authors        = refs$Authors,
  year           = suppressWarnings(as.integer(
                     sub(".*?(\\d{4}).*", "\\1", refs$Date))),
  journal        = ifelse(!is.na(refs$JournalAbbrev) & nzchar(refs$JournalAbbrev),
                          refs$JournalAbbrev, refs$JournalLong),
  volume         = as.character(refs$Vol),
  issue          = as.character(refs$Issue),
  page_start     = as.character(refs$PStart),
  page_end       = as.character(refs$PEnd),
  notes          = refs$Notes,
  stringsAsFactors = FALSE
)
for (col in c("doi","title","authors","journal","volume","issue","page_start","page_end","notes")) {
  v <- publications_df[[col]]
  publications_df[[col]][!is.na(v) & !nzchar(v)] <- NA
}
publications_df <- publications_df[publications_df$study_id %in% studies_df$study_id, ]
dbWriteTable(dst, "publications", publications_df, append = TRUE)

# --- 7. Arms --------------------------------------------------------------
# Arm universe: every (RefID, ArmNo) seen in tblArms or as a non-zero ArmNo in
# tblIntraData/tblStudyChars, restricted to kept studies.
arm_keys <- unique(rbind(
  arms_x[, c("RefID", "ArmNo")],
  intra [intra$ArmNo  > 0, c("RefID", "ArmNo")],
  chars [chars$ArmNo  > 0, c("RefID", "ArmNo")]
))
arm_keys <- arm_keys[!is.na(arm_keys$RefID) & !is.na(arm_keys$ArmNo) &
                     arm_keys$ArmNo > 0 & arm_keys$RefID %in% all_ref_ids, ]
arm_keys <- arm_keys[order(arm_keys$RefID, arm_keys$ArmNo), ]
arm_keys$arm_id <- seq_len(nrow(arm_keys))

key <- function(r, a) paste(r, a, sep = "|")

arm_name_of <- setNames(arms_x$ArmName, key(arms_x$RefID, arms_x$ArmNo))

drug_rows <- chars[chars$CatID == 57 & chars$ArmNo > 0, ]
drug_lkp  <- lookups[lookups$CatID == 57, c("AnsID", "AnsText")]
drug_name_per_arm <- setNames(
  drug_lkp$AnsText[match(drug_rows$ListVal, drug_lkp$AnsID)],
  key(drug_rows$RefID, drug_rows$ArmNo)
)

dose_rows <- chars[chars$CatID == 58 & chars$ArmNo > 0, ]
dose_amt_per_arm <- setNames(as.numeric(dose_rows$NumVal),
                             key(dose_rows$RefID, dose_rows$ArmNo))

unit_rows <- chars[chars$CatID == 59 & chars$ArmNo > 0, ]
unit_lkp  <- lookups[lookups$CatID == 59, c("AnsID", "AnsText")]
unit_name_per_arm <- setNames(
  unit_lkp$AnsText[match(unit_rows$ListVal, unit_lkp$AnsID)],
  key(unit_rows$RefID, unit_rows$ArmNo)
)

k <- key(arm_keys$RefID, arm_keys$ArmNo)
arms_df <- data.frame(
  arm_id       = arm_keys$arm_id,
  study_id     = arm_keys$RefID,
  arm_no       = arm_keys$ArmNo,
  arm_name     = unname(arm_name_of[k]),
  drug_id      = unname(drug_id_of[unname(drug_name_per_arm[k])]),
  dose_amount  = unname(dose_amt_per_arm[k]),
  dose_unit_id = unname(dose_unit_id_of[unname(unit_name_per_arm[k])]),
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "arms", arms_df, append = TRUE)

# (study_id, arm_no) -> arm_id, for the measurements join.
arm_id_of <- setNames(arms_df$arm_id, key(arms_df$study_id, arms_df$arm_no))

# --- 8. Measurements (long format) ----------------------------------------
# The 21 wide-view outcomes plus baseline PASI, main analysis only, timepoints
# within MAX_TIMEPOINT_WK week-equivalent weeks.
keep_outcome_ids <- unique(c(wide_outcomes_df$outcome_id, BASELINE_PASI_OUTCOME_ID))

m <- intra[intra$OutcomeID %in% keep_outcome_ids &
           intra$ArmNo > 0 &
           (is.na(intra$SubgroupID) | intra$SubgroupID == 0L), ]
m$subgroup_id <- 0L

# Week-equivalent timepoint cutoff. Convert the raw timepoint to weeks using the
# study's unit before comparing. Rows with no timepoint or an unknown unit are
# kept (no timepoint = not "beyond N weeks").
wk_per_unit  <- c(wk = 1, d = 1/7, mo = 52/12, min = 1/(7*24*60))
factor_for_m <- unname(wk_per_unit[unname(study_tp_unit_name[as.character(m$RefID)])])
wk_equiv     <- m$TimePeriod * factor_for_m
keep_tp <- is.na(m$TimePeriod) | is.na(factor_for_m) | wk_equiv <= MAX_TIMEPOINT_WK
m <- m[keep_tp, ]

m$arm_id <- unname(arm_id_of[key(m$RefID, m$ArmNo)])
m <- m[!is.na(m$arm_id), ]

measurements_df <- data.frame(
  arm_id      = m$arm_id,
  outcome_id  = m$OutcomeID,
  subgroup_id = m$subgroup_id,
  timepoint   = m$TimePeriod,
  k           = m$k,
  n           = m$N,
  mean        = m$Mean,
  sd          = m$SD,
  median      = m$Median,
  lo_iqr      = m$loIQR,
  hi_iqr      = m$hiIQR,
  stringsAsFactors = FALSE
)
# Collapse any accidental duplicates on (arm_id, outcome_id, subgroup_id,
# timepoint) by keeping the row with the most populated stats.
key4 <- with(measurements_df,
             paste(arm_id, outcome_id, subgroup_id, timepoint, sep = "|"))
if (anyDuplicated(key4)) {
  stat_cols <- c("k","n","mean","sd","median","lo_iqr","hi_iqr")
  score <- rowSums(!is.na(measurements_df[, stat_cols]))
  ord <- order(key4, -score)
  measurements_df <- measurements_df[ord, ]
  measurements_df <- measurements_df[!duplicated(key4[ord]), ]
}

# For each (arm_id, outcome_id, subgroup_id), keep only the largest follow-up
# timepoint. Baseline rows (timepoint == 0 or NA) are always preserved.
is_followup <- !is.na(measurements_df$timepoint) & measurements_df$timepoint > 0
fu_df  <- measurements_df[is_followup, ]
key3   <- with(fu_df, paste(arm_id, outcome_id, subgroup_id, sep = "|"))
if (anyDuplicated(key3)) {
  max_tp <- tapply(fu_df$timepoint, key3, max)
  fu_df  <- fu_df[fu_df$timepoint == max_tp[key3], ]
}
measurements_df <- rbind(measurements_df[!is_followup, ], fu_df)

dbWriteTable(dst, "measurements", measurements_df, append = TRUE)

# --- 9. Views (the contract the app reads) --------------------------------
source("R/convert/views.R")
build_views(dst)

# --- 10. Wrap up ----------------------------------------------------------
n_drug <- dbGetQuery(dst, "SELECT COUNT(DISTINCT drug) AS n FROM v_pasi WHERE drug IS NOT NULL")$n
cat(sprintf("\nDistinct drugs in v_pasi: %d\n", n_drug))
cat(sprintf("Studies: %d   Measurements: %d\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM studies")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM measurements")$n))
cat(sprintf("Publications: %d (primary: %d, secondary: %d)\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM publications")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM publications WHERE is_primary = 1")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM publications WHERE is_primary = 0")$n))

cat("\nSample v_pasi (first 5 rows):\n")
print(dbGetQuery(dst, "SELECT trial, drug, dose, timepoint, n,
                              pasi50, pasi75, pasi90, pasi100, baseline_pasi_mean
                       FROM v_pasi
                       ORDER BY trial, arm_no, timepoint
                       LIMIT 5"))

dbExecute(dst, "VACUUM")
dbDisconnect(dst)

cat(sprintf("\nDone. SQLite file: %s (%.1f MB)\n",
            sqlite_p, file.info(sqlite_p)$size / 1024 / 1024))
