# Psoriasis RCT Explorer

An R/Shiny app for browsing extracted psoriasis-RCT data (PASI, DLQI, safety
outcomes) and its meta-analysis results, over a SQLite database built from a
source Access file.

## Layout

```
.
├── app/             the Shiny app (app.R) and the sqlite db it reads
├── R/
│   ├── convert/      RevPal.accdb -> app/psoriasis-rcts.sqlite
│   ├── checks/        ad-hoc data-quality checks against the sqlite
│   └── meta-analyse/  fits the meta-analysis models, writes results into the sqlite
├── JAGS/            JAGS model definitions used by R/meta-analyse
└── claude.md        conventions for working in this repo
```

`RevPal.accdb` and `app/psoriasis-rcts.sqlite` are git-ignored — supply your
own copy of the Access database and regenerate the sqlite file locally.

## Running

All scripts assume they're run from the project root.

```powershell
Rscript R/convert/convert.R        # (re)build app/psoriasis-rcts.sqlite from RevPal.accdb
Rscript -e "shiny::runApp('app')"  # launch the app
```

The meta-analysis scripts in `R/meta-analyse/` are run manually and
separately, writing their results back into the same sqlite file.
