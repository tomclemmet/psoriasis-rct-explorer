# Psoriasis RCT Explorer

A small R/Shiny app over a normalised SQLite extract of `RevPal.accdb`. Three tabs —
**PASI**, **DLQI**, **Safety** — each with a dropdown that picks the
endpoint group shown in the table:

- **PASI**: PASI 50/75/90/100, or absolute PASI (baseline, follow-up,
  Δ from baseline on one row per arm × timepoint).
- **DLQI**: DLQI 0/1 and 0, DLQI ≤ 5, 5+/4+ point decreases, or absolute
  DLQI (baseline, follow-up, Δ).
- **Safety**: Any SAE, discontinuation (any / due to AE), serious
  infections, injection-site reactions, or malignancy (incl. NMSC).

Endpoint groups with no rows under the current filter are greyed-out in
the dropdown. Clicking a trial name opens a popover listing every
publication for that study (primary + secondaries) as Vancouver-style
citations with clickable DOIs. A summary line under each tab's dropdown
shows the current filter plus trial/publication/patient counts for the
table below. A **Download SQLite** button in the header streams the
generated `psoriasis-rcts.sqlite` file.

In the **absolute PASI** and **absolute DLQI** tables, when a study
reports baseline and follow-up but not change directly (or vice versa),
the change column is filled in as `follow-up − baseline` and marked with
`*` — a caption above the table explains. Only the mean is derived; SDs
of derived cells are left blank because they depend on the within-arm
correlation, which isn't reported.

To the left of the tables is a clickable **NMA connectivity diagram**
(built with `visNetwork`, Kamada–Kawai layout via `igraph`). Nodes are
drugs, with area proportional to the number of randomised patients
contributing to the active endpoint; edges connect drugs compared
head-to-head in the same trial, with width = number of such trials. The
network reflects the active endpoint — drugs and comparisons without
data for it are hidden. Clicking a node filters the tables to that one
drug; clicking an edge filters to the head-to-head pair. Click empty
space to clear the filter.

## Layout

```
.
├── README.md            (this file)
├── .gitignore           (excludes the .accdb and .sqlite)
├── RevPal.accdb         (not in git — supply locally; see "Cloning" below)
├── app/
│   ├── convert.R              (reads ../RevPal.accdb → writes ./psoriasis-rcts.sqlite)
│   ├── psoriasis-rcts.sqlite  (generated, not in git)
│   └── app.R            (the Shiny app — shiny::runApp("app") finds this)
└── checks/              (ad-hoc data-quality / inspection scripts; read the
                          sqlite but don't affect the app)
    ├── check.R                    (data-quality sweep over the views)
    ├── drug_doses_timepoints.R    (per-drug doses + timepoints)
    ├── publications.R             (sanity-checks the publications table)
    └── studies_multi_timepoint.R  (studies with >1 outcome timepoint)
```

`convert.R` reads the relevant `tbl*` source tables out of Access and writes
a **normalised** SQLite database: facts in `measurements` (long format, one
row per arm × outcome × subgroup × timepoint), entities in `studies` and
`arms`, and text labels in lookup tables (`drugs`, `dose_units`,
`timepoint_units`, `outcomes`, `subgroups`). It then creates four **views**
— `v_pasi`, `v_pasi_abs`, `v_dlqi`, `v_safety` — that pivot the long facts
back to the wide shape the app consumes, each keyed on
`(ref_id, arm_no, timepoint)` and restricted to `subgroup_id = 0`. The
**Download SQLite** button serves this normalised file directly; subgroup
rows (`subgroup_id > 0`) are preserved in `measurements` for downstream use.

Schema (normalised tables):

| table | role |
|---|---|
| `drugs (drug_id, drug_name)` | drug-name lookup |
| `dose_units (unit_id, unit_name)` | dose unit lookup (`mg`, `mg/kg`, …) |
| `timepoint_units (unit_id, unit_name)` | timepoint unit lookup (`wk`, `mo`) |
| `data_types (data_type_id, name)` | data-type lookup (Continuous, Dichotomous, …) from `tblDataTypes` |
| `outcomes (outcome_id, code, label, subcategory, data_type_id, endpoint_group)` | full outcome catalogue (~116 entries) from `tblOutcomeDefs`. `code` and `endpoint_group` populated only for the 21 outcomes pivoted into a view; `subcategory` distinguishes endpoint outcomes (`PASI`, `DLQI`, `Safety`) from baseline characteristics (`Demographics`, `Psoriasis characteristics`, `Previous therapy`, `Comorbidity`) |
| `subgroups (subgroup_id, subgroup_name)` | `0` = main analysis; other IDs (e.g. "Biologic-naive", "Weight < 69 kg") come from `tblSubgroups` |
| `studies (study_id, trial, design, study_date, location, phase, inclusion_criteria, exclusion_criteria, population_restriction, timepoint_unit_id)` | one row per *primary* Access RefID; `design`/`study_date`/`location`/`phase`/`inclusion_criteria`/`exclusion_criteria`/`population_restriction` come from `tblStudyChars` CatIDs 50/51/52/53/54/55/60. Secondary publications attach via `publications.study_id` rather than getting their own `studies` row |
| `publications (publication_id, study_id, is_primary, doi, title, authors, year, journal, notes)` | one row per `tblRefs` entry. `study_id` points at the primary publication's row in `studies` (the primary's own row has `study_id = publication_id` and `is_primary = 1`); secondaries link to the same primary |
| `arms (arm_id, study_id, arm_no, arm_name, drug_id, dose_amount, dose_unit_id)` | one row per study × arm; FKs into `drugs` / `dose_units` |
| `measurements (measurement_id, arm_id, outcome_id, subgroup_id, timepoint, k, n, mean, sd, median, lo_iqr, hi_iqr)` | long-format fact table; carries **every** outcome in `tblIntraData` (~4 300 rows), including baseline characteristics (age, sex, weight, BMI, ethnicity, prior therapy, baseline PASI/DLQI etc.) at `timepoint = 0`. `UNIQUE(arm_id, outcome_id, subgroup_id, timepoint)` |
| `study_subgroups (study_id, subgroup_id, n, subgroup_type, notes, excluded)` | per-study subgroup denominators from `tblSubgroupsStudies` |
| `arm_subgroups (arm_id, subgroup_id, n)` | per-arm subgroup denominators from `tblSubgroupsArms` |

All four `v_*` views share the same context columns:

| column     | source                                                       |
|------------|--------------------------------------------------------------|
| `trial`    | `studies.trial`                                              |
| `ref_id`   | `studies.study_id` (= Access RefID)                          |
| `arm_no`   | `arms.arm_no`                                                |
| `arm_name` | `arms.arm_name`                                              |
| `drug`     | `drugs.drug_name` via `arms.drug_id`                         |
| `dose`     | `arms.dose_amount` + `dose_units.unit_name` (rendered as `"40 mg"`) |
| `timepoint`| `measurements.timepoint`                                     |
| `timepoint_unit` | `timepoint_units.unit_name` via `studies.timepoint_unit_id` (usually `wk`, occasionally `mo`) |
| `n`        | `MAX(measurements.n)` across the included outcomes for that arm × timepoint |

### `v_pasi` — PASI threshold responders (binary)

| column    | OutcomeID | meaning                       |
|-----------|-----------|-------------------------------|
| `pasi50`  | 36        | `measurements.k`              |
| `pasi75`  | 13        | `measurements.k`              |
| `pasi90`  | 14        | `measurements.k`              |
| `pasi100` | 38        | `measurements.k`              |

In this schema responders are stored in column **`k`** (not `TP`) for binary outcomes.

### `v_pasi_abs` — absolute PASI score (continuous)

Pivots OutcomeID 46 (Absolute PASI) and 34 (Absolute PASI change from
baseline). Carries `_mean`, `_sd`, and `_median` for each; the app renders
`mean (SD)`.

### `v_dlqi` — DLQI endpoints (mixed)

Binary columns (k responders): `dlqi_0_1` (id 41), `dlqi_0` (51),
`dlqi_5pt_dec` (35), `dlqi_4pt_dec` (50), `dlqi_le5` (112). Continuous
columns: `abs_dlqi_*` (id 43), `abs_dlqi_change_*` (id 56).

### `v_safety` — safety outcomes (binary)

`sae` (id 20), `disc_any` (48), `disc_ae` (37), `serious_infection` (23),
`injection_site_rxn` (21), `malignancy` (96), `nmsc` (24),
`malignancy_non_nmsc` (25).

## Cloning the repo

The source Access database (`RevPal.accdb`) and the generated `psoriasis-rcts.sqlite`
are both git-ignored. After cloning, drop your own `RevPal.accdb` into the
project root and run `convert.R` (see below) to regenerate the SQLite file.

## One-time setup

R 4.5 already has `shiny`, `DBI`, `RSQLite`, `odbc`, `dplyr` installed. The app
also needs `DT`, `visNetwork`, and `meta`:

```r
install.packages(c("DT", "visNetwork", "meta"))
```

(`meta` is only needed by `meta_analyse.R` at build time. Forest plots in the
Meta-analyse modal are rendered as plain inline SVG by `app.R` itself — no
plotting library is required at runtime.)

The 64-bit "Microsoft Access Driver (*.mdb, *.accdb)" must be installed (it
already is on this machine — it ships with the 64-bit Access Database Engine).

## Run

From the project root:

```powershell
# 1. Build / rebuild the SQLite copy (only needed once, or after RevPal.accdb changes)
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" app\convert.R

# 2. Build / rebuild the meta-analysis tables (only needed once, or after
#    convert.R is re-run). Populates ma_pairwise, ma_pairwise_trials,
#    ma_proportion, ma_proportion_trials inside the same SQLite file. The
#    "Meta-analyse" button in the app reads these directly.
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" app\meta_analyse.R

# 3. Launch the Shiny app (opens in your browser)
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" -e "shiny::runApp('app', launch.browser = TRUE)"
```

## Known data quirks

- A handful of arms have `dose` populated but no `drug` (the curator filled in
  the dose without picking from the drug lookup). They show `NA` in the Drug
  column and are excluded when you filter by drug.
- A few studies have PASI outcomes recorded with no responder count (`k`); those
  cells render as blank.
- The `studies` table carries every `RefID` in the Access file (~520), but only
  the ~200 fully-extracted studies have `design` / `study_date` / inclusion-exclusion
  populated — the rest are screening-stage refs that never reached full extraction.
