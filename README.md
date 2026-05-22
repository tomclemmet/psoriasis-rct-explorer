# RevPal PASI explorer

A small R/Shiny app over a SQLite copy of `RevPal.accdb`. Shows one row per
study arm × PASI timepoint, with a multi-select **Drug** filter.

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
many `qry*` views and scratch tables) and builds a flat **`v_pasi`** table
keyed on `(ref_id, arm_no, timepoint)`.

## `v_pasi` columns

| column     | source                                                       |
|------------|--------------------------------------------------------------|
| `trial`    | `tblStudyChars.TextVal` where `CatID = 49` (Cochrane Study ID) |
| `ref_id`   | `tblIntraData.RefID`                                         |
| `arm_no`   | `tblIntraData.ArmNo`                                         |
| `arm_name` | `tblArms.ArmName`                                            |
| `drug`     | lookup of `tblStudyChars.ListVal` where `CatID = 57`         |
| `dose`     | `CatID 58` (amount, numeric) + `CatID 59` (unit, lookup)     |
| `timepoint`| `tblIntraData.TimePeriod` (weeks)                            |
| `n`        | arm denominator (max `tblIntraData.N` across the four PASI outcomes for that timepoint) |
| `pasi50`   | `tblIntraData.k` where `OutcomeID = 36`                      |
| `pasi75`   | `tblIntraData.k` where `OutcomeID = 13`                      |
| `pasi90`   | `tblIntraData.k` where `OutcomeID = 14`                      |
| `pasi100`  | `tblIntraData.k` where `OutcomeID = 38`                      |

In this schema responders are stored in column **`k`** (not `TP`) for binary outcomes.

## Cloning the repo

The source Access database (`RevPal.accdb`) and the generated `revpal.sqlite`
are both git-ignored. After cloning, drop your own `RevPal.accdb` into the
project root and run `convert.R` (see below) to regenerate the SQLite file.

## One-time setup

R 4.5 already has `shiny`, `DBI`, `RSQLite`, `odbc`, `dplyr` installed. The app
also needs `DT`:

```r
install.packages("DT")
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
