#!/usr/bin/env Rscript
# List every study whose outcome data spans more than one timepoint.
# One row per study; columns: primary_pub_id, trial, drugs, timepoints.
# Runs from Rscript, RStudio's Source button, or source() at the console.

library(DBI)
library(RSQLite)

find_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(fa)) return(dirname(normalizePath(fa[1], mustWork = FALSE)))
  for (i in seq_along(sys.frames())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  NULL
}

here <- find_script_dir()
if (is.null(here)) here <- file.path(getwd(), "checks")
sqlite_p <- normalizePath(file.path(here, "..", "app", "psoriasis-rcts.sqlite"),
                          mustWork = TRUE)

if (exists("con", envir = globalenv(), inherits = FALSE)) {
  old <- get("con", envir = globalenv())
  if (inherits(old, "DBIConnection") && DBI::dbIsValid(old)) {
    try(DBI::dbDisconnect(old), silent = TRUE)
  }
  rm("con", envir = globalenv())
}

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqlite_p)
stopifnot(DBI::dbIsValid(con))

# Pull one row per (study, drug, timepoint) seen in the main-analysis
# outcome data, then aggregate in R.
rows <- DBI::dbGetQuery(con, "
  SELECT DISTINCT
         s.study_id   AS primary_pub_id,
         s.trial      AS trial,
         dr.drug_name AS drug,
         m.timepoint  AS timepoint,
         tu.unit_name AS timepoint_unit
  FROM   measurements m
  JOIN   arms a              ON a.arm_id = m.arm_id
  JOIN   studies s           ON s.study_id = a.study_id
  LEFT   JOIN drugs dr       ON dr.drug_id = a.drug_id
  LEFT   JOIN timepoint_units tu ON tu.unit_id = s.timepoint_unit_id
  WHERE  m.subgroup_id = 0
    AND  m.timepoint IS NOT NULL
    AND  m.timepoint <> 0
    AND  m.timepoint <= 24
")

DBI::dbDisconnect(con)

# Keep studies with >1 distinct timepoint.
tp_counts <- tapply(rows$timepoint, rows$primary_pub_id,
                    function(x) length(unique(x)))
multi_ids <- as.integer(names(tp_counts)[tp_counts > 1])
rows <- rows[rows$primary_pub_id %in% multi_ids, ]

# Build one row per study.
out <- do.call(rbind, lapply(sort(unique(rows$primary_pub_id)), function(id) {
  sub  <- rows[rows$primary_pub_id == id, ]
  unit <- unique(sub$timepoint_unit[!is.na(sub$timepoint_unit)])
  unit <- if (length(unit) == 1) unit else ""
  tps  <- sort(unique(sub$timepoint))
  data.frame(
    primary_pub_id = id,
    trial          = unique(sub$trial)[1],
    drugs          = paste(sort(unique(sub$drug[!is.na(sub$drug)])),
                           collapse = ", "),
    timepoints     = paste0(paste(tps, collapse = ", "),
                            if (nzchar(unit)) paste0(" ", unit) else ""),
    stringsAsFactors = FALSE
  )
}))

cat(sprintf("%d studies have outcome data at >1 timepoint.\n\n", nrow(out)))
print(out, row.names = FALSE, right = FALSE)
