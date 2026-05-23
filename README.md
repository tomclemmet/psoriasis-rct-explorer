# Psoriasis RCT Explorer

A small R/Shiny app over a SQLite copy of `RevPal.accdb`. Three tabs —
**PASI**, **DLQI**, **Safety** — each with a dropdown that picks the
endpoint group shown in the table:

- **PASI**: PASI 50/75/90/100, or absolute PASI (baseline, follow-up,
  Δ from baseline on one row per arm × timepoint).
- **DLQI**: DLQI 0/1 and 0, DLQI ≤ 5, 5+/4+ point decreases, or absolute
  DLQI (baseline, follow-up, Δ).
- **Safety**: Any SAE, discontinuation (any / due to AE), serious
  infections, injection-site reactions, or malignancy (incl. NMSC).

Endpoint groups with no rows under the current filter are greyed-out in
the dropdown. Trial cells link to the publication DOI when one is
recorded in `tblRefs` (about 12% of refs in this dataset have a DOI).
A **Download SQLite** button in the filter bar streams the generated
`revpal.sqlite` file.

To the left of the tables is a clickable **NMA connectivity diagram**
(built with `visNetwork`, Kamada–Kawai layout via `igraph`). Nodes are
drugs, sized by the number of trials in which they appear; edges connect
drugs compared head-to-head in the same trial, with width = number of
such trials. Clicking a node filters the tables to that one drug;
clicking an edge filters to the head-to-head pair (only trials including
both drugs). Click empty space, or the **Clear** button, to reset.

## Layout

```
.
├── README.md            (this file)
├── .gitignore           (excludes the .accdb and .sqlite)
├── RevPal.accdb         (not in git — supply locally; see "Cloning" below)
└── app/
    ├── convert.R        (reads ../RevPal.accdb → writes ./revpal.sqlite)
    ├── revpal.sqlite    (generated, not in git)
    └── app.R            (the Shiny app — shiny::runApp("app") finds this)
```

`convert.R` copies the real `tbl*` source tables out of Access (dropping the
many `qry*` views and scratch tables) and builds four flat **`v_*`** tables,
each keyed on `(ref_id, arm_no, timepoint)` and restricted to
`SubgroupID = 0` (main arms only).

All four views share the same context columns:

| column     | source                                                       |
|------------|--------------------------------------------------------------|
| `trial`    | `tblStudyChars.TextVal` where `CatID = 49` (Cochrane Study ID) |
| `ref_id`   | `tblIntraData.RefID`                                         |
| `arm_no`   | `tblIntraData.ArmNo`                                         |
| `arm_name` | `tblArms.ArmName`                                            |
| `drug`     | lookup of `tblStudyChars.ListVal` where `CatID = 57`         |
| `dose`     | `CatID 58` (amount, numeric) + `CatID 59` (unit, lookup)     |
| `timepoint`| `tblIntraData.TimePeriod` (numeric)                          |
| `timepoint_unit` | `tblLongitudinalDataDefs.Unit` → `tblLongitudinalUnitDefs.strUnit` (usually `wk`, occasionally `mo`) |
| `n`        | arm denominator (max `tblIntraData.N` across the included outcomes for that timepoint) |

### `v_pasi` — PASI threshold responders (binary)

| column    | OutcomeID | meaning                       |
|-----------|-----------|-------------------------------|
| `pasi50`  | 36        | `tblIntraData.k`              |
| `pasi75`  | 13        | `tblIntraData.k`              |
| `pasi90`  | 14        | `tblIntraData.k`              |
| `pasi100` | 38        | `tblIntraData.k`              |

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

The source Access database (`RevPal.accdb`) and the generated `revpal.sqlite`
are both git-ignored. After cloning, drop your own `RevPal.accdb` into the
project root and run `convert.R` (see below) to regenerate the SQLite file.

## One-time setup

R 4.5 already has `shiny`, `DBI`, `RSQLite`, `odbc`, `dplyr` installed. The app
also needs `DT` and `visNetwork`:

```r
install.packages(c("DT", "visNetwork"))
```

The 64-bit "Microsoft Access Driver (*.mdb, *.accdb)" must be installed (it
already is on this machine — it ships with the 64-bit Access Database Engine).

## Run

From the project root:

```powershell
# 1. Build / rebuild the SQLite copy (only needed once, or after RevPal.accdb changes)
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" app\convert.R

# 2. Launch the Shiny app (opens in your browser)
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" -e "shiny::runApp('app', launch.browser = TRUE)"
```

## Known data quirks

- A handful of arms have `dose` populated but no `drug` (the curator filled in
  the dose without picking from the drug lookup). They show `NA` in the Drug
  column and are excluded when you filter by drug.
- A few studies have PASI outcomes recorded with no responder count (`k`); those
  cells render as blank.
