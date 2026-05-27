#!/usr/bin/env Rscript
# Convert RevPal.accdb -> psoriasis-rcts.sqlite as a normalised relational database,
# then expose four wide SQL views (v_pasi, v_pasi_abs, v_dlqi, v_safety) that
# the Shiny app reads. Run from the project root:
#   Rscript app/convert.R
#
# Output schema (one fact per row, no repeated text):
#
#   Lookups
#     drugs(drug_id PK, drug_name UNIQUE)
#     dose_units(unit_id PK, unit_name UNIQUE)
#     timepoint_units(unit_id PK, unit_name UNIQUE)
#     data_types(data_type_id PK, name)            -- from tblDataTypes
#     outcomes(outcome_id PK, code, label, subcategory,
#              data_type_id FK, endpoint_group)    -- full Access catalogue
#                                                  -- (incl. baseline chars)
#     subgroups(subgroup_id PK, subgroup_name)
#
#   Entities
#     studies(study_id PK, trial,
#             design, study_date, location, phase,
#             inclusion_criteria, exclusion_criteria,
#             population_restriction,
#             timepoint_unit_id FK)
#         -- Keyed by the *primary* publication's RefID. Secondary
#         -- publications are recorded in `publications` and never appear
#         -- here (the source data attributes all measurements/arms/chars
#         -- to the primary RefID).
#     publications(publication_id PK, study_id FK, is_primary,
#                  doi, title, authors, year, journal, notes)
#         -- One row per tblRefs entry; secondaries link to their primary
#         -- via `study_id` (= tblStudyDefs.ParentID). The primary's row
#         -- has study_id == publication_id and is_primary = 1.
#     arms(arm_id PK, study_id FK, arm_no, arm_name,
#          drug_id FK, dose_amount, dose_unit_id FK)
#
#   Facts
#     measurements(measurement_id PK, arm_id FK, outcome_id FK,
#                  subgroup_id FK, timepoint,
#                  k, n, mean, sd, median, lo_iqr, hi_iqr)
#         -- carries ALL outcomes from tblIntraData, not just the four
#         -- endpoint groups. Demographics / baseline characteristics
#         -- (age, sex, weight, BMI, baseline PASI, ethnicity, prior
#         -- therapy etc.) live here as their own OutcomeIDs.
#     study_subgroups(study_id, subgroup_id, n,
#                     subgroup_type, notes, excluded)
#     arm_subgroups(arm_id, subgroup_id, n)
#
#   Views (rebuild the wide tables the app expects, filtered to main arms)
#     v_pasi, v_pasi_abs, v_dlqi, v_safety
#
# Subgroup rows are preserved in `measurements` (subgroup_id > 0); the four
# views filter to subgroup_id = 0 so the app behaves as before.

suppressPackageStartupMessages({
  library(DBI)
  library(odbc)
  library(RSQLite)
})

resolve_script_dir <- function() {
  # Rscript: --file=... is on commandArgs.
  file_arg <- sub("^--file=", "",
                  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) && nzchar(file_arg[[1]])) {
    return(normalizePath(dirname(file_arg[[1]]), mustWork = TRUE))
  }
  # RStudio: ask the IDE for the active source editor path.
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (nzchar(p)) return(normalizePath(dirname(p), mustWork = TRUE))
  }
  # source() from the console: sys.frames() carries the ofile.
  for (i in rev(seq_along(sys.frames()))) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(normalizePath(dirname(ofile), mustWork = TRUE))
  }
  # Fallback: assume cwd is the project root.
  normalizePath("app", mustWork = TRUE)
}

here     <- resolve_script_dir()
accdb    <- normalizePath(file.path(here, "..", "RevPal.accdb"), mustWork = TRUE)
sqlite_p <- file.path(here, "psoriasis-rcts.sqlite")

# Close any lingering connection to the target from a previous run in this
# R session - otherwise file.remove() silently fails on Windows and the
# subsequent dbConnect() can hand back an invalid handle.
if (exists("dst", envir = globalenv(), inherits = FALSE)) {
  try(dbDisconnect(get("dst", envir = globalenv())), silent = TRUE)
  rm("dst", envir = globalenv())
}
if (exists("src", envir = globalenv(), inherits = FALSE)) {
  try(dbDisconnect(get("src", envir = globalenv())), silent = TRUE)
  rm("src", envir = globalenv())
}

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
# names the view. Only the 21 outcomes used by a view appear here; the full
# 116-entry catalogue (which adds demographics / baseline characteristics etc.)
# is read from tblOutcomeDefs below and merged with this decoration.
wide_outcomes_df <- read.table(text = "
outcome_id|code|endpoint_group
36|pasi50|pasi
13|pasi75|pasi
14|pasi90|pasi
38|pasi100|pasi
46|abs_pasi|pasi_abs
34|abs_pasi_change|pasi_abs
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
sgnames <- dbReadTable(src, "tblSubgroups")  # SubgroupID, SubgroupName, intSubGOrder
odefs   <- dbReadTable(src, "tblOutcomeDefs")    # 116 outcomes catalogue
dtypes  <- dbReadTable(src, "tblDataTypes")      # DataTypeID -> name
sg_st   <- dbReadTable(src, "tblSubgroupsStudies") # per (study, subgroup)
sg_arm  <- dbReadTable(src, "tblSubgroupsArms")    # per (study, subgroup, arm)
studefs <- dbReadTable(src, "tblStudyDefs")        # RefID, ParentID, Notes, ...
dbDisconnect(src)

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
CREATE TABLE timepoint_units (
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
  endpoint_group  TEXT                   -- 'pasi','pasi_abs','dlqi','safety',
                                         -- or NULL when not exposed by a view
);
CREATE TABLE subgroups (
  subgroup_id    INTEGER PRIMARY KEY,
  subgroup_name  TEXT
);
CREATE TABLE studies (
  study_id                INTEGER PRIMARY KEY,  -- = primary publication's RefID
  trial                   TEXT,    -- Cochrane Study ID (CatID 49)
  design                  TEXT,    -- CatID 50
  study_date              TEXT,    -- CatID 51 (free-form date range)
  location                TEXT,    -- CatID 52
  phase                   TEXT,    -- CatID 53 (resolved from lookup)
  inclusion_criteria      TEXT,    -- CatID 54
  exclusion_criteria      TEXT,    -- CatID 55
  population_restriction  TEXT,    -- CatID 60
  timepoint_unit_id       INTEGER REFERENCES timepoint_units(unit_id)
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
  volume          TEXT,                        -- tblRefs.Vol (text: occasionally non-numeric)
  issue           TEXT,                        -- tblRefs.Issue
  page_start      TEXT,                        -- tblRefs.PStart (text: e.g. S12)
  page_end        TEXT,                        -- tblRefs.PEnd
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
CREATE TABLE study_subgroups (
  study_id       INTEGER NOT NULL REFERENCES studies(study_id),
  subgroup_id    INTEGER NOT NULL REFERENCES subgroups(subgroup_id),
  n              INTEGER,
  subgroup_type  TEXT,
  notes          TEXT,
  excluded       INTEGER,  -- 0/1 from tblSubgroupsStudies.blnExc
  PRIMARY KEY (study_id, subgroup_id)
);
CREATE TABLE arm_subgroups (
  arm_id         INTEGER NOT NULL REFERENCES arms(arm_id),
  subgroup_id    INTEGER NOT NULL REFERENCES subgroups(subgroup_id),
  n              INTEGER,
  PRIMARY KEY (arm_id, subgroup_id)
);
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

# Full outcome catalogue from tblOutcomeDefs, decorated with view metadata
# for the 21 outcomes the app pivots into v_pasi / v_pasi_abs / v_dlqi /
# v_safety.
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

# Subgroup IDs referenced anywhere: measurements, per-study denominators,
# per-arm denominators. Plus 0 (the main-analysis sentinel) which has no
# row in tblSubgroups but is implied throughout.
subgroup_ids <- sort(unique(as.integer(c(
  0L, intra$SubgroupID, sg_st$SubgroupID, sg_arm$SubgroupID
))))
subgroup_ids <- subgroup_ids[!is.na(subgroup_ids)]
sg_name_of <- setNames(sgnames$SubgroupName, as.character(sgnames$SubgroupID))
subgroups_df <- data.frame(
  subgroup_id   = subgroup_ids,
  subgroup_name = ifelse(subgroup_ids == 0L,
                         "Main analysis",
                         unname(sg_name_of[as.character(subgroup_ids)])),
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "subgroups", subgroups_df, append = TRUE)

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

# Timepoint units. Include "wk" as the default even if absent in source.
tp_unit_names <- sort(unique(c("wk", lunits$strUnit)))
tp_unit_names <- tp_unit_names[!is.na(tp_unit_names) & nzchar(tp_unit_names)]
timepoint_units_df <- data.frame(
  unit_id   = seq_along(tp_unit_names),
  unit_name = tp_unit_names,
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "timepoint_units", timepoint_units_df, append = TRUE)
tp_unit_id_of <- setNames(timepoint_units_df$unit_id, timepoint_units_df$unit_name)

# --- 6. Studies -----------------------------------------------------------
# Universe = *primary* RefIDs that have any extracted data. tblStudyDefs
# distinguishes primaries (ParentID NULL/0) from secondaries (ParentID > 0);
# every measurement/arm/char row in the source attributes data to a primary
# RefID, so secondaries never need a `studies` row. Their bibliography info
# lives in `publications` further down.
primary_ids <- studefs$RefID[is.na(studefs$ParentID) | studefs$ParentID == 0]

# Fail loud if a secondary has measurements or arms attributed directly —
# would silently drop data once we restrict to primaries below.
# (tblStudyChars also records CatID 49 = trial name redundantly against
# secondaries; that's a known denormalisation in the source and is ignored.)
secondary_with_data <- setdiff(
  unique(as.integer(c(intra$RefID, arms_x$RefID))),
  primary_ids
)
secondary_with_data <- secondary_with_data[!is.na(secondary_with_data)]
if (length(secondary_with_data)) {
  stop("Secondary RefIDs have measurements or arms: ",
       paste(sort(secondary_with_data), collapse = ", "))
}

# Universe: primaries that have any genuine data. tblStudyChars rows
# against secondaries (CatID 49 only — redundant trial name) are excluded.
data_ref_ids <- unique(as.integer(c(
  intra$RefID,
  chars$RefID[chars$RefID %in% primary_ids],
  arms_x$RefID
)))
data_ref_ids <- data_ref_ids[!is.na(data_ref_ids)]
all_ref_ids <- sort(intersect(primary_ids, data_ref_ids))

# Helper: pull a study-level TextVal field (ArmNo = 0) for a given CatID
# into a named vector keyed by RefID.
study_text_of <- function(cat_id) {
  r <- chars[chars$CatID == cat_id &
             (is.na(chars$ArmNo) | chars$ArmNo == 0), ]
  setNames(r$TextVal, as.character(r$RefID))
}

# Study Phase (CatID 53) is a ListVal — resolve via tblStudyCharsDefsLookups.
phase_rows <- chars[chars$CatID == 53 &
                    (is.na(chars$ArmNo) | chars$ArmNo == 0), ]
phase_lkp  <- lookups[lookups$CatID == 53, c("AnsID", "AnsText")]
phase_of   <- setNames(
  phase_lkp$AnsText[match(phase_rows$ListVal, phase_lkp$AnsID)],
  as.character(phase_rows$RefID)
)

trial_of      <- study_text_of(49)
design_of     <- study_text_of(50)
date_of       <- study_text_of(51)
location_of   <- study_text_of(52)
inc_of        <- study_text_of(54)
exc_of        <- study_text_of(55)
pop_restr_of  <- study_text_of(60)

# Per-study timepoint unit. Use MIN(strUnit) across every outcome present
# in tblLongitudinalDataDefs for the study; deterministic, constant per
# study in this dataset.
ld_join <- merge(
  lddefs[, c("RefID", "Unit")],
  lunits[, c("UnitID", "strUnit")],
  by.x = "Unit", by.y = "UnitID", all.x = TRUE
)
tp_unit_by_ref <- tapply(ld_join$strUnit, ld_join$RefID,
                         function(x) min(x[!is.na(x)]))

pick <- function(map, ids) unname(map[as.character(ids)])

studies_df <- data.frame(
  study_id               = all_ref_ids,
  trial                  = pick(trial_of,      all_ref_ids),
  design                 = pick(design_of,     all_ref_ids),
  study_date             = pick(date_of,       all_ref_ids),
  location               = pick(location_of,   all_ref_ids),
  phase                  = pick(phase_of,      all_ref_ids),
  inclusion_criteria     = pick(inc_of,        all_ref_ids),
  exclusion_criteria     = pick(exc_of,        all_ref_ids),
  population_restriction = pick(pop_restr_of,  all_ref_ids),
  timepoint_unit_id      = unname(tp_unit_id_of[
    ifelse(is.na(tp_unit_by_ref[as.character(all_ref_ids)]),
           "wk", tp_unit_by_ref[as.character(all_ref_ids)])
  ]),
  stringsAsFactors = FALSE
)
# Blank strings -> NA.
for (col in c("trial","design","study_date","location","phase",
              "inclusion_criteria","exclusion_criteria","population_restriction")) {
  v <- studies_df[[col]]
  studies_df[[col]][!is.na(v) & !nzchar(v)] <- NA
}
dbWriteTable(dst, "studies", studies_df, append = TRUE)

# --- 6b. Publications -----------------------------------------------------
# One row per tblRefs entry. Each publication points at its primary's
# study_id (= ParentID if non-NULL/>0, else its own RefID). Only keep
# publications whose resolved primary actually has a `studies` row, so the
# FK is satisfied; secondaries dangling off an excluded primary are dropped
# with a warning.
parent_of <- setNames(
  ifelse(is.na(studefs$ParentID) | studefs$ParentID == 0,
         studefs$RefID, studefs$ParentID),
  as.character(studefs$RefID)
)
pub_study_id <- unname(parent_of[as.character(refs$ID)])
# Publications missing from tblStudyDefs default to "primary of themselves".
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
# Blank strings -> NA in text columns.
for (col in c("doi","title","authors","journal","volume","issue","page_start","page_end","notes")) {
  v <- publications_df[[col]]
  publications_df[[col]][!is.na(v) & !nzchar(v)] <- NA
}

dropped_pubs <- sum(!(publications_df$study_id %in% studies_df$study_id))
if (dropped_pubs) {
  cat(sprintf("  publications: dropping %d row(s) whose primary is outside the study universe\n",
              dropped_pubs))
}
publications_df <- publications_df[publications_df$study_id %in% studies_df$study_id, ]
dbWriteTable(dst, "publications", publications_df, append = TRUE)

# --- 7. Arms --------------------------------------------------------------
# Arm universe: every (RefID, ArmNo) seen in tblArms or as a non-zero ArmNo
# in tblIntraData/tblStudyChars.
arm_keys <- unique(rbind(
  arms_x[, c("RefID", "ArmNo")],
  intra [intra$ArmNo  > 0, c("RefID", "ArmNo")],
  chars [chars$ArmNo  > 0, c("RefID", "ArmNo")]
))
arm_keys <- arm_keys[!is.na(arm_keys$RefID) & !is.na(arm_keys$ArmNo) &
                     arm_keys$ArmNo > 0, ]
arm_keys <- arm_keys[order(arm_keys$RefID, arm_keys$ArmNo), ]
arm_keys$arm_id <- seq_len(nrow(arm_keys))

# Lookups by (RefID, ArmNo) key.
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
drug_names_resolved <- unname(drug_name_per_arm[k])
unit_names_resolved <- unname(unit_name_per_arm[k])

arms_df <- data.frame(
  arm_id       = arm_keys$arm_id,
  study_id     = arm_keys$RefID,
  arm_no       = arm_keys$ArmNo,
  arm_name     = unname(arm_name_of[k]),
  drug_id      = unname(drug_id_of[drug_names_resolved]),
  dose_amount  = unname(dose_amt_per_arm[k]),
  dose_unit_id = unname(dose_unit_id_of[unit_names_resolved]),
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "arms", arms_df, append = TRUE)

# (arm_id, study_id, arm_no) -> arm_id, for the measurements join.
arm_id_of <- setNames(arms_df$arm_id, key(arms_df$study_id, arms_df$arm_no))

# --- 8. Measurements (long format, all outcomes, all subgroups) -----------
# Carry every tblIntraData row to measurements. Demographics / baseline
# characteristics (Age, Sex, Weight, BMI, baseline PASI, ethnicity, prior
# therapy, etc.) land here alongside the on-treatment endpoints. The four
# v_* views still filter to the 21 wide-view outcomes, so the app behaves
# identically. The %in% guard is just FK safety — drops rows whose
# OutcomeID is missing from tblOutcomeDefs (none in this dataset).
m <- intra[intra$OutcomeID %in% outcomes_df$outcome_id & intra$ArmNo > 0, ]
m$subgroup_id <- ifelse(is.na(m$SubgroupID), 0L, as.integer(m$SubgroupID))
m$arm_id <- unname(arm_id_of[key(m$RefID, m$ArmNo)])

if (any(is.na(m$arm_id))) {
  bad <- m[is.na(m$arm_id), c("RefID","ArmNo","OutcomeID")][1:5, , drop = FALSE]
  warning("Some intra rows did not match an arm: ", nrow(m[is.na(m$arm_id),]))
  print(bad)
  m <- m[!is.na(m$arm_id), ]
}

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
# timepoint) by taking the row with the most populated stats.
key4 <- with(measurements_df,
             paste(arm_id, outcome_id, subgroup_id, timepoint, sep = "|"))
if (anyDuplicated(key4)) {
  stat_cols <- c("k","n","mean","sd","median","lo_iqr","hi_iqr")
  score <- rowSums(!is.na(measurements_df[, stat_cols]))
  ord <- order(key4, -score)
  measurements_df <- measurements_df[ord, ]
  measurements_df <- measurements_df[!duplicated(
    key4[ord]
  ), ]
}
dbWriteTable(dst, "measurements", measurements_df, append = TRUE)

# --- 8b. Subgroup denominators --------------------------------------------
study_sg <- sg_st[sg_st$RefID %in% all_ref_ids &
                  sg_st$SubgroupID %in% subgroup_ids, ]
study_subgroups_df <- data.frame(
  study_id      = study_sg$RefID,
  subgroup_id   = study_sg$SubgroupID,
  n             = as.integer(study_sg$N),
  subgroup_type = study_sg$SubgroupType,
  notes         = study_sg$SubgroupNotes,
  excluded      = as.integer(as.logical(study_sg$blnExc)),
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "study_subgroups", study_subgroups_df, append = TRUE)

arm_sg <- sg_arm[sg_arm$SubgroupID %in% subgroup_ids, ]
arm_sg$arm_id <- unname(arm_id_of[key(arm_sg$RefID, arm_sg$ArmNo)])
arm_sg <- arm_sg[!is.na(arm_sg$arm_id), ]
arm_subgroups_df <- data.frame(
  arm_id      = arm_sg$arm_id,
  subgroup_id = arm_sg$SubgroupID,
  n           = as.integer(arm_sg$N),
  stringsAsFactors = FALSE
)
dbWriteTable(dst, "arm_subgroups", arm_subgroups_df, append = TRUE)

# --- 9. Views (the contract the app reads) --------------------------------
# Common arm/study context CTE used by every view. Mirrors the column names
# the app already expects: trial, ref_id, arm_no, arm_name, drug, dose,
# timepoint_unit. Dose is rebuilt as "<amount> <unit>" (integer-valued
# amounts are rendered without a trailing .0, preserving the prior format).
arm_ctx_cte <- "
WITH arm_ctx AS (
  SELECT a.arm_id,
         s.study_id     AS ref_id,
         s.trial        AS trial,
         a.arm_no       AS arm_no,
         a.arm_name     AS arm_name,
         dr.drug_name   AS drug,
         TRIM(
           CASE
             WHEN a.dose_amount IS NULL THEN ''
             WHEN a.dose_amount = CAST(a.dose_amount AS INTEGER)
               THEN CAST(CAST(a.dose_amount AS INTEGER) AS TEXT)
             ELSE CAST(a.dose_amount AS TEXT)
           END
           || ' ' || COALESCE(du.unit_name, '')
         )              AS dose,
         tu.unit_name   AS timepoint_unit
  FROM   arms a
  JOIN   studies s         ON s.study_id = a.study_id
  LEFT   JOIN drugs dr     ON dr.drug_id = a.drug_id
  LEFT   JOIN dose_units du ON du.unit_id = a.dose_unit_id
  LEFT   JOIN timepoint_units tu ON tu.unit_id = s.timepoint_unit_id
)"

# Build one view. `pivot_cols` is a vector of "MAX(CASE WHEN o.code = 'X' ...) AS col"
# clauses; `outcome_codes` lists the codes the view consults (drives the N
# computation and the measurements filter).
build_view <- function(name, outcome_codes, pivot_cols, select_cols) {
  cat(sprintf("Building view %s ...\n", name))
  in_list <- paste(sprintf("'%s'", outcome_codes), collapse = ", ")
  pivots  <- paste(pivot_cols, collapse = ",\n           ")
  selects <- paste(select_cols, collapse = ", ")

  sql <- sprintf("
    CREATE VIEW %s AS
    %s,
    pivot AS (
      SELECT m.arm_id, m.timepoint,
             %s,
             MAX(m.n) AS n
      FROM   measurements m
      JOIN   outcomes o ON o.outcome_id = m.outcome_id
      WHERE  o.code IN (%s)
        AND  m.subgroup_id = 0
      GROUP BY m.arm_id, m.timepoint
    )
    SELECT ctx.trial,
           ctx.ref_id,
           ctx.arm_no,
           ctx.arm_name,
           ctx.drug,
           ctx.dose,
           pivot.timepoint AS timepoint,
           ctx.timepoint_unit,
           pivot.n         AS n,
           %s
    FROM   pivot
    JOIN   arm_ctx ctx ON ctx.arm_id = pivot.arm_id
  ", name, arm_ctx_cte, pivots, in_list, selects)

  dbExecute(dst, sql)
  n_rows <- dbGetQuery(dst, sprintf("SELECT COUNT(*) AS n FROM %s", name))$n
  cat(sprintf("  %s rows: %d\n", name, n_rows))
}

# Helper to build a MAX(CASE WHEN ... THEN <field> END) AS <alias> clause.
case_k    <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.k      END) AS %s", code, alias)
case_mean <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.mean   END) AS %s", code, alias)
case_sd   <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.sd     END) AS %s", code, alias)
case_med  <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.median END) AS %s", code, alias)
case_lo   <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.lo_iqr END) AS %s", code, alias)
case_hi   <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.hi_iqr END) AS %s", code, alias)

# v_pasi - PASI threshold responders.
build_view(
  "v_pasi",
  outcome_codes = c("pasi50","pasi75","pasi90","pasi100"),
  pivot_cols = c(case_k("pasi50",  "pasi50"),
                 case_k("pasi75",  "pasi75"),
                 case_k("pasi90",  "pasi90"),
                 case_k("pasi100", "pasi100")),
  select_cols = c("pivot.pasi50", "pivot.pasi75", "pivot.pasi90", "pivot.pasi100")
)

# v_pasi_abs - absolute / change-from-baseline PASI.
build_view(
  "v_pasi_abs",
  outcome_codes = c("abs_pasi","abs_pasi_change"),
  pivot_cols = c(case_mean("abs_pasi",         "abs_pasi_mean"),
                 case_sd  ("abs_pasi",         "abs_pasi_sd"),
                 case_med ("abs_pasi",         "abs_pasi_median"),
                 case_lo  ("abs_pasi",         "abs_pasi_lo_iqr"),
                 case_hi  ("abs_pasi",         "abs_pasi_hi_iqr"),
                 case_mean("abs_pasi_change",  "abs_pasi_change_mean"),
                 case_sd  ("abs_pasi_change",  "abs_pasi_change_sd"),
                 case_med ("abs_pasi_change",  "abs_pasi_change_median")),
  select_cols = c("pivot.abs_pasi_mean","pivot.abs_pasi_sd","pivot.abs_pasi_median",
                  "pivot.abs_pasi_lo_iqr","pivot.abs_pasi_hi_iqr",
                  "pivot.abs_pasi_change_mean","pivot.abs_pasi_change_sd",
                  "pivot.abs_pasi_change_median")
)

# v_dlqi - DLQI binary + continuous endpoints.
build_view(
  "v_dlqi",
  outcome_codes = c("dlqi_0_1","dlqi_0","dlqi_5pt_dec","dlqi_4pt_dec","dlqi_le5",
                    "abs_dlqi","abs_dlqi_change"),
  pivot_cols = c(case_k   ("dlqi_0_1",        "dlqi_0_1"),
                 case_k   ("dlqi_0",          "dlqi_0"),
                 case_k   ("dlqi_5pt_dec",    "dlqi_5pt_dec"),
                 case_k   ("dlqi_4pt_dec",    "dlqi_4pt_dec"),
                 case_k   ("dlqi_le5",        "dlqi_le5"),
                 case_mean("abs_dlqi",        "abs_dlqi_mean"),
                 case_sd  ("abs_dlqi",        "abs_dlqi_sd"),
                 case_med ("abs_dlqi",        "abs_dlqi_median"),
                 case_mean("abs_dlqi_change", "abs_dlqi_change_mean"),
                 case_sd  ("abs_dlqi_change", "abs_dlqi_change_sd"),
                 case_med ("abs_dlqi_change", "abs_dlqi_change_median")),
  select_cols = c("pivot.dlqi_0_1","pivot.dlqi_0","pivot.dlqi_5pt_dec",
                  "pivot.dlqi_4pt_dec","pivot.dlqi_le5",
                  "pivot.abs_dlqi_mean","pivot.abs_dlqi_sd","pivot.abs_dlqi_median",
                  "pivot.abs_dlqi_change_mean","pivot.abs_dlqi_change_sd",
                  "pivot.abs_dlqi_change_median")
)

# v_safety - binary safety outcomes.
build_view(
  "v_safety",
  outcome_codes = c("sae","disc_any","disc_ae","serious_infection",
                    "injection_site_rxn","malignancy","nmsc","malignancy_non_nmsc"),
  pivot_cols = c(case_k("sae",                 "sae"),
                 case_k("disc_any",            "disc_any"),
                 case_k("disc_ae",             "disc_ae"),
                 case_k("serious_infection",   "serious_infection"),
                 case_k("injection_site_rxn",  "injection_site_rxn"),
                 case_k("malignancy",          "malignancy"),
                 case_k("nmsc",                "nmsc"),
                 case_k("malignancy_non_nmsc", "malignancy_non_nmsc")),
  select_cols = c("pivot.sae","pivot.disc_any","pivot.disc_ae",
                  "pivot.serious_infection","pivot.injection_site_rxn",
                  "pivot.malignancy","pivot.nmsc","pivot.malignancy_non_nmsc")
)

# --- 10. Wrap up ----------------------------------------------------------
n_drug <- dbGetQuery(dst, "SELECT COUNT(DISTINCT drug) AS n FROM v_pasi WHERE drug IS NOT NULL")$n
cat(sprintf("\nDistinct drugs in v_pasi: %d\n", n_drug))
cat(sprintf("Outcomes catalogued: %d  (of which exposed via a view: %d)\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM outcomes")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM outcomes WHERE code IS NOT NULL")$n))
cat(sprintf("Measurements rows (all outcomes, all subgroups): %d\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM measurements")$n))
cat(sprintf("  of which subgroup_id > 0: %d\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM measurements WHERE subgroup_id > 0")$n))
cat(sprintf("  baseline-characteristic rows (Demographics / Psoriasis chars / Previous therapy / Comorbidity): %d\n",
            dbGetQuery(dst, "
              SELECT COUNT(*) AS n FROM measurements m
              JOIN outcomes o ON o.outcome_id = m.outcome_id
              WHERE o.subcategory IN ('Demographics','Psoriasis characteristics',
                                      'Previous therapy','Comorbidity')")$n))
cat(sprintf("study_subgroups rows: %d   arm_subgroups rows: %d\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM study_subgroups")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM arm_subgroups")$n))
cat(sprintf("Studies with non-null design: %d / %d\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM studies WHERE design IS NOT NULL")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM studies")$n))
cat(sprintf("Publications: %d (primary: %d, secondary: %d)\n",
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM publications")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM publications WHERE is_primary = 1")$n,
            dbGetQuery(dst, "SELECT COUNT(*) AS n FROM publications WHERE is_primary = 0")$n))

cat("\nSample v_pasi (first 5 rows):\n")
print(dbGetQuery(dst, "SELECT trial, drug, dose, timepoint, n, pasi50, pasi75, pasi90, pasi100
                       FROM v_pasi
                       ORDER BY trial, arm_no, timepoint
                       LIMIT 5"))

dbExecute(dst, "VACUUM")
dbDisconnect(dst)

cat(sprintf("\nDone. SQLite file: %s (%.1f MB)\n",
            sqlite_p, file.info(sqlite_p)$size / 1024 / 1024))
