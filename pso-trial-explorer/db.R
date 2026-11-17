# One read-only connection shared by the whole app process (the db is never
# written to at runtime), closed when the app stops.
DB_CON <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
onStop(function() dbDisconnect(DB_CON))

read_db <- function(sql, params = list()) {
  if (length(params)) dbGetQuery(DB_CON, sql, params = params) else dbGetQuery(DB_CON, sql)
}

DB_TABLES <- dbListTables(DB_CON)

# The per-arm views behind the data tables are small (a few hundred rows
# each), so they're read once at startup and filtered in memory rather than
# re-queried on every filter change. Sorted like the old SQL
# "ORDER BY trial, arm_no, timepoint" (byte-order text, NULLs first).
order_like_sql <- function(df, cols) {
  keys <- unname(as.list(df[cols]))
  df <- df[do.call(order, c(keys, list(na.last = FALSE, method = "radix"))), , drop = FALSE]
  rownames(df) <- NULL
  df
}

VIEW_DATA <- lapply(
  c(v_pasi = "v_pasi", v_dlqi = "v_dlqi", v_safety = "v_safety"),
  function(v) order_like_sql(read_db(sprintf("SELECT * FROM %s", v)),
                             c("trial", "arm_no", "timepoint"))
)

query_view <- function(table, state) {
  df <- VIEW_DATA[[table]]
  keep <- if (is.null(state)) {
    rep(TRUE, nrow(df))
  } else if (identical(state$kind, "node")) {
    df$drug %in% state$drug
  } else if (identical(state$kind, "edge")) {
    both <- intersect(df$trial[df$drug %in% state$from],
                      df$trial[df$drug %in% state$to])
    df$drug %in% c(state$from, state$to) & df$trial %in% both
  } else {
    rep(TRUE, nrow(df))
  }
  df <- df[keep, , drop = FALSE]
  rownames(df) <- NULL
  df
}

# Whether the meta-analysis view exists yet. It's built by
# R/analysis/app-meta-analysis.R, a separate (slow, JAGS-driven) step from
# R/convert/trial-estimates.R's v_trial_estimates -- right after reconverting
# the Access db, trial estimates are often present while this isn't.
# fetch_ma() and fetch_ma_axis_combos() then return 0-row results so forest
# plots still render from trial estimates alone, with no diamonds.
HAS_MA     <- "v_meta_analysis"   %in% DB_TABLES
HAS_TRIALS <- "v_trial_estimates" %in% DB_TABLES

# Meta-analysis results and per-trial estimates are read into memory once
# and split by endpoint, so building a forest plot is a handful of data-frame
# subsets rather than (previously) up to ~150 SQL round trips per modal.
MA_ALL <- if (HAS_MA) read_db("SELECT * FROM v_meta_analysis") else NULL
MA_BY_ENDPOINT <- if (HAS_MA) split(MA_ALL, MA_ALL$endpoint) else list()

TRIALS_BY_ENDPOINT <- if (HAS_TRIALS) {
  te <- read_db("SELECT * FROM v_trial_estimates")
  split(te, te$endpoint)
} else list()

EMPTY_MA_RESULT <- data.frame(
  endpoint = character(), type = character(), effects = character(),
  comp_tx = character(), ref_tx = character(), measure = character(),
  method = character(), likelihood = character(),
  mean = numeric(), lower = numeric(), upper = numeric(), dic = numeric(),
  stringsAsFactors = FALSE
)

# Subset `df` to rows where each named column in `filters` matches its
# value(s). NULL filters are skipped; a single NA matches missing (SQL "IS
# NULL"). %in% rather than == so NA cells never match a non-NA value.
filter_rows <- function(df, filters) {
  keep <- rep(TRUE, nrow(df))
  for (col in names(filters)) {
    v <- filters[[col]]
    if (is.null(v)) next
    keep <- keep & if (length(v) == 1 && is.na(v)) is.na(df[[col]]) else df[[col]] %in% v
  }
  out <- df[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Columns treated as "which model produced this row" axes, dynamically
# discoverable and filterable rather than hardcoded — see MA_AXIS_COLUMNS in
# app.R and the project memory: meta-analysis-modelling-seam.
MA_DISTINCT_COLUMNS <- c("method", "likelihood", "metareg")

# `axes` is a named list filtering "which model produced this row" columns
# (see MA_DISTINCT_COLUMNS) — e.g. list(likelihood = "multinomial",
# method = "class_effects").
fetch_ma <- function(endpoint, type = NULL, effects = NULL,
                     comp_tx = NULL, ref_tx = NULL, measure = NULL, axes = list()) {
  df <- MA_BY_ENDPOINT[[endpoint]]
  if (is.null(df)) return(EMPTY_MA_RESULT)
  stopifnot(all(names(axes) %in% MA_DISTINCT_COLUMNS))
  filter_rows(df, c(list(type = type, effects = effects, comp_tx = comp_tx,
                         ref_tx = ref_tx, measure = measure), axes))
}

fetch_trials <- function(endpoint, comp_tx = NULL, ref_tx = NULL, measure = NULL) {
  df <- TRIALS_BY_ENDPOINT[[endpoint]]
  if (is.null(df)) return(data.frame())
  filter_rows(df, list(comp_tx = comp_tx, ref_tx = ref_tx, measure = measure))
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

# Distinct combinations of several axis columns present for a set of
# endpoints + a result type (one row per combination actually occurring in
# the data, not the full cartesian product of each column's individual
# values) — drives dynamic multi-axis discovery in app.R.
fetch_ma_axis_combos <- function(endpoint, type, columns, effects = NULL) {
  stopifnot(all(columns %in% MA_DISTINCT_COLUMNS))
  if (!HAS_MA) return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
  df <- filter_rows(MA_ALL[c("endpoint", "type", "effects", columns)],
                    list(endpoint = endpoint, type = type, effects = effects))
  out <- unique(df[, columns, drop = FALSE])
  rownames(out) <- NULL
  out
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
  lapply(c(pasi = "v_pasi", dlqi = "v_dlqi", safety = "v_safety"), function(v) {
    df <- VIEW_DATA[[v]]
    order_like_sql(df[df$ref_id %in% study_id, , drop = FALSE],
                   c("arm_no", "timepoint"))
  })
}
