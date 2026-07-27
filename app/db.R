read_db <- function(sql, params = list()) {
  con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  if (length(params)) dbGetQuery(con, sql, params = params) else dbGetQuery(con, sql)
}

query_view <- function(table, state) {
  base_order <- "ORDER BY trial, arm_no, timepoint"
  if (is.null(state)) {
    return(read_db(sprintf("SELECT * FROM %s %s", table, base_order)))
  }
  if (identical(state$kind, "node")) {
    return(read_db(
      sprintf("SELECT * FROM %s WHERE drug = ? %s", table, base_order),
      params = list(state$drug)
    ))
  }
  if (identical(state$kind, "edge")) {
    return(read_db(
      sprintf("SELECT * FROM %s
               WHERE drug IN (?, ?)
                 AND trial IN (
                   SELECT trial FROM %s WHERE drug = ?
                   INTERSECT
                   SELECT trial FROM %s WHERE drug = ?
                 )
               %s", table, table, table, base_order),
      params = list(state$from, state$to, state$from, state$to)
    ))
  }
  read_db(sprintf("SELECT * FROM %s %s", table, base_order))
}

# `axes` is a named list filtering "which model produced this row" columns
# (see MA_DISTINCT_COLUMNS below) — e.g. list(likelihood = "multinomial",
# method = "class_effects"). Column names are checked against that allowlist
# since they're interpolated into the SQL (not bind parameters).
fetch_ma <- function(endpoint, type = NULL, effects = NULL,
                     comp_tx = NULL, ref_tx = NULL, measure = NULL, axes = list()) {
  conds  <- "WHERE endpoint = ?"
  params <- list(endpoint)
  if (!is.null(type))    { conds <- paste(conds, "AND type = ?");    params <- c(params, list(type)) }
  if (!is.null(effects)) { conds <- paste(conds, "AND effects = ?"); params <- c(params, list(effects)) }
  if (!is.null(comp_tx)) { conds <- paste(conds, "AND comp_tx = ?"); params <- c(params, list(comp_tx)) }
  if (!is.null(ref_tx))  { conds <- paste(conds, "AND ref_tx = ?");  params <- c(params, list(ref_tx)) }
  if (!is.null(measure)) { conds <- paste(conds, "AND measure = ?"); params <- c(params, list(measure)) }
  for (col in names(axes)) {
    stopifnot(col %in% MA_DISTINCT_COLUMNS)
    conds  <- paste(conds, sprintf("AND %s = ?", col))
    params <- c(params, list(axes[[col]]))
  }
  read_db(sprintf("SELECT * FROM v_meta_analysis %s", conds), params = params)
}

fetch_trials <- function(endpoint, comp_tx = NULL, ref_tx = NULL, measure = NULL) {
  conds  <- "WHERE endpoint = ?"
  params <- list(endpoint)
  if (!is.null(comp_tx)) { conds <- paste(conds, "AND comp_tx = ?"); params <- c(params, list(comp_tx)) }
  if (!is.null(ref_tx))  { conds <- paste(conds, "AND ref_tx = ?");  params <- c(params, list(ref_tx)) }
  if (!is.null(measure)) { conds <- paste(conds, "AND measure = ?"); params <- c(params, list(measure)) }
  read_db(sprintf("SELECT * FROM v_trial_estimates %s", conds), params = params)
}

fetch_ma_directed <- function(endpoint, type, effects, comp, ref, measure, axes = list()) {
  r <- fetch_ma(endpoint, type = type, effects = effects,
                comp_tx = comp, ref_tx = ref, measure = measure, axes = axes)
  if (nrow(r)) return(list(mean = r$mean[1], lower = r$lower[1], upper = r$upper[1], dic = r$dic[1]))
  r <- fetch_ma(endpoint, type = type, effects = effects,
                comp_tx = ref, ref_tx = comp, measure = measure, axes = axes)
  if (nrow(r)) return(list(mean = -r$mean[1], lower = -r$upper[1], upper = -r$lower[1], dic = r$dic[1]))
  NULL
}

# Columns treated as "which model produced this row" axes, dynamically
# discoverable and filterable rather than hardcoded — see MA_AXIS_COLUMNS in
# app.R and the project memory: meta-analysis-modelling-seam.
MA_DISTINCT_COLUMNS <- c("method", "likelihood")

# Distinct combinations of several axis columns present for a set of
# endpoints + a result type (one row per combination actually occurring in
# the data, not the full cartesian product of each column's individual
# values) — drives dynamic multi-axis discovery in app.R.
fetch_ma_axis_combos <- function(endpoint, type, columns) {
  stopifnot(all(columns %in% MA_DISTINCT_COLUMNS))
  in_list  <- paste(rep("?", length(endpoint)), collapse = ", ")
  cols_sql <- paste(columns, collapse = ", ")
  sql <- sprintf(
    "SELECT DISTINCT %s FROM v_meta_analysis WHERE type = ? AND endpoint IN (%s)",
    cols_sql, in_list
  )
  read_db(sql, params = c(list(type), as.list(endpoint)))
}

coalesce0 <- function(x) ifelse(is.na(x), 0L, as.integer(x))

fetch_baselines <- function(study_id) {
  read_db(
    "SELECT o.label, o.subcategory, a.arm_no, a.arm_name,
            d.drug_name AS drug,
            TRIM(
              TRIM(
                CASE
                  WHEN a.dose_amount IS NULL THEN ''
                  WHEN a.dose_amount = CAST(a.dose_amount AS INTEGER)
                    THEN CAST(CAST(a.dose_amount AS INTEGER) AS TEXT)
                  ELSE CAST(a.dose_amount AS TEXT)
                END
                || ' ' || COALESCE(du.unit_name, '')
              )
              || ' ' || COALESCE(fr.frequency_name, '')
            ) AS dose,
            m.n, m.k, m.mean, m.sd
     FROM   measurements m
     JOIN   outcomes o      ON o.outcome_id = m.outcome_id
     JOIN   arms a          ON a.arm_id     = m.arm_id
     LEFT   JOIN drugs d    ON d.drug_id    = a.drug_id
     LEFT   JOIN dose_units du ON du.unit_id = a.dose_unit_id
     LEFT   JOIN frequencies fr ON fr.frequency_id = a.frequency_id
     WHERE  a.study_id = ?
       AND  o.code IS NULL
       AND  o.subcategory IN ('Demographics', 'Psoriasis characteristics',
                              'Previous therapy', 'Comorbidity')
     ORDER  BY o.subcategory, o.outcome_id, a.arm_no",
    params = list(study_id)
  )
}

fetch_trial_results <- function(study_id) {
  list(
    pasi   = read_db("SELECT * FROM v_pasi   WHERE ref_id = ? ORDER BY arm_no, timepoint",
                     params = list(study_id)),
    dlqi   = read_db("SELECT * FROM v_dlqi   WHERE ref_id = ? ORDER BY arm_no, timepoint",
                     params = list(study_id)),
    safety = read_db("SELECT * FROM v_safety WHERE ref_id = ? ORDER BY arm_no, timepoint",
                     params = list(study_id))
  )
}
