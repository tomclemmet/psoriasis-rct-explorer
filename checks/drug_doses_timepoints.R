#!/usr/bin/env Rscript
# For each drug in psoriasis-rcts.sqlite, print the distinct doses and
# timepoints observed across all arms / measurements. Runs from Rscript,
# RStudio's Source button, or source() at the console.

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

# Close anything stale from a previous run in this R session.
if (exists("con", envir = globalenv(), inherits = FALSE)) {
  try(DBI::dbDisconnect(get("con", envir = globalenv())), silent = TRUE)
  rm("con", envir = globalenv())
}

cat("RSQLite:", as.character(utils::packageVersion("RSQLite")),
    "  DBI:",  as.character(utils::packageVersion("DBI")),
    "\nFile:  ", sqlite_p, "\n")

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqlite_p)
stopifnot(DBI::dbIsValid(con))

rows <- DBI::dbGetQuery(con, "
  SELECT dr.drug_name                            AS drug,
         a.dose_amount                           AS dose_amount,
         du.unit_name                            AS dose_unit,
         m.timepoint                             AS timepoint,
         tu.unit_name                            AS timepoint_unit
  FROM   arms a
  JOIN   drugs dr            ON dr.drug_id = a.drug_id
  LEFT   JOIN dose_units du  ON du.unit_id = a.dose_unit_id
  JOIN   studies s           ON s.study_id = a.study_id
  LEFT   JOIN timepoint_units tu ON tu.unit_id = s.timepoint_unit_id
  LEFT   JOIN measurements m ON m.arm_id  = a.arm_id
                            AND m.subgroup_id = 0
  ORDER BY dr.drug_name, a.dose_amount, m.timepoint
")

DBI::dbDisconnect(con)

fmt_dose <- function(amount, unit) {
  if (is.na(amount)) return(NA_character_)
  txt <- if (amount == as.integer(amount)) format(as.integer(amount))
         else format(amount)
  if (is.na(unit) || !nzchar(unit)) txt else paste(txt, unit)
}

for (drug in sort(unique(rows$drug))) {
  sub <- rows[rows$drug == drug, ]
  doses <- unique(mapply(fmt_dose, sub$dose_amount, sub$dose_unit,
                         USE.NAMES = FALSE))
  doses <- doses[!is.na(doses)]
  tps   <- sort(unique(sub$timepoint[!is.na(sub$timepoint)]))
  tp_unit <- unique(sub$timepoint_unit[!is.na(sub$timepoint_unit)])
  tp_unit <- if (length(tp_unit) == 1) tp_unit else ""

  cat(drug, "\n", sep = "")
  cat("  doses:      ",
      if (length(doses)) paste(doses, collapse = ", ") else "(none recorded)",
      "\n", sep = "")
  cat("  timepoints: ",
      if (length(tps)) paste0(paste(tps, collapse = ", "),
                              if (nzchar(tp_unit)) paste0(" ", tp_unit) else "")
      else "(none recorded)",
      "\n\n", sep = "")
}
