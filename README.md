# Psoriasis RCT Explorer

Repository for a network meta-analysis project in psoriasis, including R and JAGS code used to analyse data and an interactive R shiny app to present the results, available at [shiny.its.manchester.ac.uk/pso-trial-explorer/](https://shiny.its.manchester.ac.uk/pso-trial-explorer/). Trial results and meta-analysis estimates are stored in a sqlite database that is queried by the Shiny app. For further details see the about section of the Shiny app.

## Layout

```
.
├───pso-trial-explorer/            The R shiny app containing data extraction from trials and meta-analysis results
├───JAGS/                          JAGS code used in the multinomial network meta-analyses
└───R/
    ├───analysis/                  Meta-analysis code
    │   ├───bivariate/             Bivariate network meta-analysis code (in development)
    │   └───metaregression/        Meta-regression analyses
    ├───convert/                   Code to process raw data into the correct format for the Shiny app
    ├───data-checks/               Checks for the data extraction
    └───utils/                     Helper functions for data processing, plots, and model fitting
```


