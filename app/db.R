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

fetch_ma <- function(endpoint, type = NULL, effects = NULL,
                     comp_tx = NULL, ref_tx = NULL, measure = NULL, method = NULL) {
  conds  <- "WHERE endpoint = ?"
  params <- list(endpoint)
  if (!is.null(type))    { conds <- paste(conds, "AND type = ?");    params <- c(params, list(type)) }
  if (!is.null(effects)) { conds <- paste(conds, "AND effects = ?"); params <- c(params, list(effects)) }
  if (!is.null(comp_tx)) { conds <- paste(conds, "AND comp_tx = ?"); params <- c(params, list(comp_tx)) }
  if (!is.null(ref_tx))  { conds <- paste(conds, "AND ref_tx = ?");  params <- c(params, list(ref_tx)) }
  if (!is.null(measure)) { conds <- paste(conds, "AND measure = ?"); params <- c(params, list(measure)) }
  if (!is.null(method))  { conds <- paste(conds, "AND method = ?");  params <- c(params, list(method)) }
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

fetch_ma_directed <- function(endpoint, type, effects, comp, ref, measure, method = NULL) {
  r <- fetch_ma(endpoint, type = type, effects = effects,
                comp_tx = comp, ref_tx = ref, measure = measure, method = method)
  if (nrow(r)) return(list(mean = r$mean[1], lower = r$lower[1], upper = r$upper[1]))
  r <- fetch_ma(endpoint, type = type, effects = effects,
                comp_tx = ref, ref_tx = comp, measure = measure, method = method)
  if (nrow(r)) return(list(mean = -r$mean[1], lower = -r$upper[1], upper = -r$lower[1]))
  NULL
}

coalesce0 <- function(x) ifelse(is.na(x), 0L, as.integer(x))

fetch_baselines <- function(study_id) {
  read_db(
    "SELECT o.label, o.subcategory, a.arm_no, a.arm_name,
            d.drug_name AS drug,
            TRIM(
              CASE
                WHEN a.dose_amount IS NULL THEN ''
                WHEN a.dose_amount = CAST(a.dose_amount AS INTEGER)
                  THEN CAST(CAST(a.dose_amount AS INTEGER) AS TEXT)
                ELSE CAST(a.dose_amount AS TEXT)
              END
              || ' ' || COALESCE(du.unit_name, '')
            ) AS dose,
            m.n, m.k, m.mean, m.sd
     FROM   measurements m
     JOIN   outcomes o      ON o.outcome_id = m.outcome_id
     JOIN   arms a          ON a.arm_id     = m.arm_id
     LEFT   JOIN drugs d    ON d.drug_id    = a.drug_id
     LEFT   JOIN dose_units du ON du.unit_id = a.dose_unit_id
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
