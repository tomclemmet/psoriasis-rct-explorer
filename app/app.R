suppressPackageStartupMessages({
  library(shiny)
  library(DBI)
  library(RSQLite)
  library(DT)
  library(visNetwork)
  library(jsonlite)
  library(ggplot2)
  library(ggiraph)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

DB_PATH <- file.path(dirname(sys.frame(1)$ofile %||% "."), "psoriasis-rcts.sqlite")
if (!file.exists(DB_PATH)) DB_PATH <- "app/psoriasis-rcts.sqlite"
if (!file.exists(DB_PATH)) stop("psoriasis-rcts.sqlite not found - run convert.R first.")

read_db <- function(sql, params = list()) {
  con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  if (length(params)) dbGetQuery(con, sql, params = params) else dbGetQuery(con, sql)
}

# Turn a per-arm (trial, drug, ref_id, arm_no, n_arm) data frame into
# visNetwork node + edge data frames. Nodes = drugs, sized by total patients
# across surviving arms (area-proportional via sqrt). Edges = unordered drug
# pairs that appear together in a trial, width = number of trials with that
# head-to-head. Pure function: called once at startup to build the master
# graph (union across all endpoint groups) and again per-group for subsets.
NODE_SIZE_MIN <- 14
NODE_SIZE_MAX <- 55

# Scale a vector of patient counts to a visible-radius range using sqrt so
# that area ∝ patients. `ref_max` anchors the upper end so the same patient
# count produces the same on-screen size across every endpoint subset —
# without this, visNetwork's per-render auto-scaling would make 1,000
# patients look the same as 100 in views where 100 is the largest.
size_from_patients <- function(n, ref_max) {
  if (!length(n)) return(numeric(0))
  if (!is.finite(ref_max) || ref_max <= 0) return(rep(NODE_SIZE_MIN, length(n)))
  r <- sqrt(pmax(n, 1)) / sqrt(ref_max)
  pmin(NODE_SIZE_MAX, pmax(NODE_SIZE_MIN,
                            NODE_SIZE_MIN + (NODE_SIZE_MAX - NODE_SIZE_MIN) * r))
}

build_network_data <- function(td, ref_max_n = NA_real_) {
  td <- td[!is.na(td$drug) & nzchar(td$drug) & !is.na(td$trial), , drop = FALSE]
  # One row per (trial, arm) keeps the patient count from being multi-counted
  # if the same arm produced several rows in the source view.
  td <- unique(td[, c("trial", "ref_id", "arm_no", "drug", "n_arm"),
                  drop = FALSE])

  empty_edges <- data.frame(from = character(), to = character(),
                            n_trials = integer(), id = character(),
                            value = integer(), title = character(),
                            stringsAsFactors = FALSE)
  if (!nrow(td)) {
    return(list(
      nodes = data.frame(id = character(), label = character(),
                         value = numeric(), title = character(),
                         stringsAsFactors = FALSE),
      edges = empty_edges
    ))
  }

  drugs <- sort(unique(td$drug))
  n_patients <- vapply(drugs, function(d)
    sum(td$n_arm[td$drug == d], na.rm = TRUE), numeric(1))
  n_trials_per_drug <- vapply(drugs, function(d)
    length(unique(td$trial[td$drug == d])), integer(1))

  size_ref <- if (is.finite(ref_max_n) && ref_max_n > 0) ref_max_n
              else max(n_patients, 1, na.rm = TRUE)
  nodes_df <- data.frame(
    id    = drugs,
    label = drugs,
    # Absolute size in pixels (radius). Bypassing visNetwork's `scaling.value`
    # auto-rescale so a given patient count looks identical across endpoint
    # subsets. sqrt → area ∝ patients (standard NMA-diagram convention).
    size  = size_from_patients(n_patients, size_ref),
    title = sprintf("<b>%s</b><br/>%d trial(s), %s patient(s)",
                    drugs, n_trials_per_drug,
                    formatC(n_patients, format = "d", big.mark = ",")),
    stringsAsFactors = FALSE
  )
  nodes_df <- nodes_df[order(-n_patients[match(nodes_df$id, drugs)],
                             nodes_df$id), ]

  td_pair <- unique(td[, c("trial", "drug"), drop = FALSE])
  pair_rows <- list()
  for (tr in unique(td_pair$trial)) {
    ds <- sort(unique(td_pair$drug[td_pair$trial == tr]))
    if (length(ds) < 2) next
    cmb <- utils::combn(ds, 2)
    pair_rows[[tr]] <- data.frame(from = cmb[1, ], to = cmb[2, ],
                                  stringsAsFactors = FALSE)
  }
  if (!length(pair_rows)) {
    return(list(nodes = nodes_df, edges = empty_edges))
  }
  pairs <- do.call(rbind, pair_rows)
  pair_counts <- as.data.frame(table(pairs$from, pairs$to),
                               stringsAsFactors = FALSE)
  names(pair_counts) <- c("from", "to", "n_trials")
  pair_counts <- pair_counts[pair_counts$n_trials > 0, ]
  # Stable id derived from (from, to) so the same comparison keeps the same
  # edge id across subset rebuilds — that lets visNetworkProxy reconcile
  # additions/removals correctly when the user switches endpoints.
  pair_counts$id    <- sprintf("e_%s__%s", pair_counts$from, pair_counts$to)
  pair_counts$value <- pair_counts$n_trials
  pair_counts$title <- sprintf("<b>%s &harr; %s</b><br/>%d trial(s)",
                               pair_counts$from, pair_counts$to,
                               pair_counts$n_trials)

  list(nodes = nodes_df, edges = pair_counts)
}

# study_id -> list of publication rows (primary first, then secondaries by
# year). Built once at startup; used by fmt_citations() to render the
# popover that opens when a user clicks a trial name.
.pub_rows <- read_db(
  "SELECT study_id, is_primary, doi, title, authors, year, journal,
          volume, issue, page_start, page_end
   FROM   publications
   ORDER  BY study_id, is_primary DESC, year ASC"
)
pubs_by_study <- split(.pub_rows, .pub_rows$study_id)

# Pre-render one citations-HTML string per study_id at startup, then
# attribute-escape it for embedding in data-citations="...". A trial
# typically appears in dozens of rows (one per arm × timepoint), so doing
# this per-row in fmt_trial() was the hot path slowing table renders.
# Built lazily below once fmt_citations() is defined.

# Reformat one author entry from "Griffiths,C.E." to "Griffiths CE".
.fmt_one_author <- function(a) {
  a <- trimws(a)
  if (!nzchar(a)) return(NA_character_)
  # First comma separates surname from initials.
  parts <- strsplit(a, ",", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(a)
  surname  <- trimws(parts[1])
  initials <- paste(parts[-1], collapse = ",")
  initials <- gsub("[.[:space:]]+", "", initials)
  paste(surname, initials)
}

# Reformat a semicolon-separated author string to Vancouver style:
# "Griffiths,C.E.; Strober,B.E." -> "Griffiths CE, Strober BE".
# Truncates to first 6 authors + ", et al" when >6 (standard Vancouver).
.fmt_authors <- function(s) {
  if (is.na(s) || !nzchar(s)) return("")
  authors <- vapply(strsplit(s, ";", fixed = TRUE)[[1]],
                    .fmt_one_author, character(1))
  authors <- authors[!is.na(authors) & nzchar(authors)]
  if (!length(authors)) return("")
  if (length(authors) > 6) {
    paste0(paste(authors[1:6], collapse = ", "), ", et al")
  } else {
    paste(authors, collapse = ", ")
  }
}

# Build the "Year;Vol(Iss):Pstart-Pend" fragment. Any missing piece drops
# out cleanly. Returns "" if year is also missing (caller decides).
.fmt_volissue <- function(year, volume, issue, p_start, p_end) {
  out <- if (is.na(year)) "" else as.character(year)
  if (!is.na(volume) && nzchar(volume)) {
    out <- paste0(out, ";", volume)
    if (!is.na(issue) && nzchar(issue)) out <- paste0(out, "(", issue, ")")
  }
  if (!is.na(p_start) && nzchar(p_start)) {
    out <- paste0(out, ":", p_start)
    if (!is.na(p_end) && nzchar(p_end) && p_end != p_start) {
      out <- paste0(out, "-", p_end)
    }
  }
  out
}

# Render one publications data frame as an HTML <ol> of Vancouver-style
# citations. Returns "" when no publications are known for the study.
fmt_citations <- function(pubs) {
  esc <- htmltools::htmlEscape
  if (is.null(pubs) || !nrow(pubs)) return("")
  items <- character(nrow(pubs))
  for (i in seq_len(nrow(pubs))) {
    p <- pubs[i, ]
    authors <- .fmt_authors(p$authors)
    title   <- if (is.na(p$title)) "" else trimws(p$title)
    journal <- if (is.na(p$journal)) "" else trimws(p$journal)
    # Some "publications" are registry entries with a URL in the journal
    # field; render those as a link and skip vol/iss/pages.
    is_url <- nzchar(journal) && grepl("^https?://", journal)
    parts <- c()
    if (nzchar(authors)) parts <- c(parts, paste0(esc(authors), "."))
    if (nzchar(title))   parts <- c(parts, paste0(esc(title), "."))
    if (is_url) {
      parts <- c(parts, sprintf(
        '<a href="%s" target="_blank" rel="noopener">%s</a>.',
        esc(journal), esc(journal)))
    } else {
      if (nzchar(journal)) parts <- c(parts, paste0(esc(journal), "."))
      vi <- .fmt_volissue(p$year, p$volume, p$issue, p$page_start, p$page_end)
      if (nzchar(vi)) parts <- c(parts, paste0(esc(vi), "."))
    }
    if (!is.na(p$doi) && nzchar(p$doi)) {
      parts <- c(parts, sprintf(
        '<a href="https://doi.org/%s" target="_blank" rel="noopener">https://doi.org/%s</a>',
        esc(p$doi), esc(p$doi)))
    }
    items[i] <- paste0('<div class="trial-cite">', paste(parts, collapse = " "), "</div>")
  }
  paste0('<div class="trial-cites">', paste(items, collapse = ""), "</div>")
}

# Pre-rendered citation HTML per study_id. Held as plain (unescaped) HTML
# because it's delivered to the browser once as a JS object (window.studyCitations)
# rather than being embedded into a `data-` attribute on every table cell.
# fmt_trial() only writes the study_id into each cell; the JS click handler
# looks the HTML up at click time.
cites_html_by_study <- vapply(
  pubs_by_study, fmt_citations, character(1)
)
cites_html_by_study <- cites_html_by_study[nzchar(cites_html_by_study)]
studies_with_cites <- names(cites_html_by_study)
cites_json <- jsonlite::toJSON(as.list(cites_html_by_study), auto_unbox = TRUE)

# Baseline PASI recorded as a "Psoriasis characteristics" outcome
# (outcome_id 11) — used as a fallback for the Absolute PASI table when an
# arm has no week-0 abs_pasi row. Curators usually record baseline PASI as
# a baseline characteristic; week 0 of the absolute-PASI longitudinal series
# is only present for the minority of arms where that timepoint was
# explicitly extracted.
.baseline_pasi <- read_db("
  SELECT a.study_id AS ref_id, a.arm_no AS arm_no,
         MAX(m.mean) AS mean, MAX(m.sd) AS sd
  FROM   measurements m
  JOIN   arms a ON a.arm_id = m.arm_id
  WHERE  m.outcome_id = 11 AND m.subgroup_id = 0
  GROUP  BY a.study_id, a.arm_no
")
.baseline_pasi_key <- paste(.baseline_pasi$ref_id, .baseline_pasi$arm_no, sep = "|")
baseline_pasi_lookup <- list(
  mean = setNames(.baseline_pasi$mean, .baseline_pasi_key),
  sd   = setNames(.baseline_pasi$sd,   .baseline_pasi_key)
)

# Render trial text as a popover trigger. Clicking the trial name opens a
# small floating panel listing every publication for that study (primary +
# secondaries) as Vancouver-style citations with clickable DOIs. The popover
# is implemented in vanilla JS (see tags$script in the UI) via event
# delegation, so it survives DT redraws without a drawCallback.
#
# Each anchor carries only `data-ref-id`; the JS handler looks the citation
# HTML up in window.studyCitations. Trials with no publications fall back to
# a span (no popover).
# Caller must pass `escape = FALSE` for the Trial column.
fmt_trial <- function(trial, ref_id) {
  labels <- htmltools::htmlEscape(trial)
  ids    <- as.character(ref_id)
  has    <- ids %in% studies_with_cites
  out    <- paste0("<span>", labels, "</span>")
  out[has] <- sprintf(
    '<a href="javascript:void(0)" class="trial-pop" data-ref-id="%s">%s</a>',
    ids[has], labels[has]
  )
  out
}

# Filter state shape:
#   NULL                                        - no filter, show all rows
#   list(kind = "node", drug = "Adalimumab")    - single-drug filter
#   list(kind = "edge", from = "A", to = "B")   - head-to-head pair: only
#                                                 rows where drug in {A,B}
#                                                 AND trial includes both
#                                                 drugs (per v_pasi).
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
                   SELECT trial FROM v_pasi WHERE drug = ?
                   INTERSECT
                   SELECT trial FROM v_pasi WHERE drug = ?
                 )
               %s", table, base_order),
      params = list(state$from, state$to, state$from, state$to)
    ))
  }
  read_db(sprintf("SELECT * FROM %s %s", table, base_order))
}

fmt_pasi <- function(k, n) {
  out <- rep("", length(k))
  ok  <- !is.na(k) & !is.na(n) & n > 0
  pct <- round(k[ok] / n[ok] * 100)
  out[ok] <- sprintf("%d (%d%%)", as.integer(k[ok]), pct)
  out
}

fmt_mean_sd <- function(mean, sd, digits = 1) {
  out <- rep("", length(mean))
  ok  <- !is.na(mean)
  m   <- formatC(mean[ok], format = "f", digits = digits)
  s   <- ifelse(is.na(sd[ok]), "",
                sprintf(" (%s)", formatC(sd[ok], format = "f", digits = digits)))
  out[ok] <- paste0(m, s)
  out
}

# Fill in the change column when both baseline and follow-up are reported
# but the study didn't report change directly. Returns the (possibly
# updated) change vector and a logical flag per row marking derived cells.
# Only the mean is derived: the SD of the difference depends on the
# within-arm baseline/follow-up correlation, which is rarely reported, so
# SDs of derived cells are left blank rather than guessed.
derive_change <- function(baseline, follow, change) {
  cd <- rep(FALSE, length(change))
  m <- is.na(change) & !is.na(baseline) & !is.na(follow)
  change[m] <- follow[m] - baseline[m]
  cd[m] <- TRUE
  list(change = change, change_derived = cd)
}

# fmt_mean_sd, but derived cells are rendered without an SD (we don't have
# one) and tagged with a trailing asterisk so the table caption can explain
# them. Non-derived cells render exactly as fmt_mean_sd would.
fmt_mean_sd_marked <- function(mean, sd, derived, digits = 1) {
  sd_eff <- sd
  sd_eff[derived] <- NA_real_
  out <- fmt_mean_sd(mean, sd_eff, digits = digits)
  mark <- derived & nzchar(out)
  out[mark] <- paste0(out[mark], "*")
  out
}

fmt_timepoint <- function(timepoint, unit) {
  unit_lbl <- ifelse(unit == "wk", "wks", unit)
  ifelse(is.na(timepoint), "", paste(timepoint, unit_lbl))
}

# Build the Drug cell text: "Adalimumab 40 mg, 16 wks". Dose/timepoint
# omitted when missing.
fmt_drug <- function(drug, dose, timepoint, unit) {
  tp_txt   <- fmt_timepoint(timepoint, unit)
  has_dose <- !is.na(dose) & nzchar(dose)
  with_dose <- ifelse(has_dose, paste0(drug, " ", dose), drug)
  ifelse(nzchar(tp_txt), paste0(with_dose, ", ", tp_txt), with_dose)
}

# Per-arm baseline: value at timepoint == 0 for the same (ref_id, arm_no),
# with optional fallback to a pre-loaded `list(mean=<named>, sd=<named>)`
# keyed by "ref_id|arm_no" when no timepoint-0 row exists. Returns vectors
# aligned with df rows.
baseline_lookup <- function(df, mean_col, sd_col, fallback = NULL) {
  key <- paste(df$ref_id, df$arm_no, sep = "|")
  is_b <- !is.na(df$timepoint) & df$timepoint == 0
  bkey <- paste(df$ref_id[is_b], df$arm_no[is_b], sep = "|")
  i    <- match(key, bkey)
  mean <- df[[mean_col]][is_b][i]
  sd   <- df[[sd_col]][is_b][i]
  if (!is.null(fallback)) {
    miss <- is.na(mean)
    mean[miss] <- unname(fallback$mean[key[miss]])
    sd[miss]   <- unname(fallback$sd[key[miss]])
  }
  list(mean = mean, sd = sd)
}

# Generic binary-subset formatter: format each column as "k (pct%)", build
# the Drug cell, drop rows with no data in any of the selected endpoints,
# return Trial/Drug/N + selected columns.
format_binary_subset <- function(df, cols) {
  for (col in cols) df[[col]] <- fmt_pasi(df[[col]], df$n)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  has_any <- Reduce(`|`, lapply(cols, function(c) nzchar(df[[c]])))
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", cols)]
}

format_pasi_response <- function(df) {
  format_binary_subset(df, c("pasi50", "pasi75", "pasi90", "pasi100"))
}

format_pasi_absolute <- function(df) {
  b <- baseline_lookup(df, "abs_pasi_mean", "abs_pasi_sd",
                       fallback = baseline_pasi_lookup)
  d <- derive_change(b$mean, df$abs_pasi_mean, df$abs_pasi_change_mean)
  df$baseline        <- fmt_mean_sd(b$mean, b$sd)
  df$on_tx           <- fmt_mean_sd(df$abs_pasi_mean, df$abs_pasi_sd)
  df$abs_pasi_change <- fmt_mean_sd_marked(d$change, df$abs_pasi_change_sd,
                                           d$change_derived)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  # Drop pure baseline rows; each follow-up row now carries its baseline.
  df <- df[is.na(df$timepoint) | df$timepoint > 0, ]
  has_any <- nzchar(df$baseline) | nzchar(df$on_tx) | nzchar(df$abs_pasi_change)
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", "baseline", "on_tx", "abs_pasi_change")]
}

format_dlqi_zero <- function(df) {
  format_binary_subset(df, c("dlqi_0_1", "dlqi_0"))
}

format_dlqi_threshold <- function(df) {
  format_binary_subset(df, c("dlqi_le5"))
}

format_dlqi_change <- function(df) {
  format_binary_subset(df, c("dlqi_5pt_dec", "dlqi_4pt_dec"))
}

format_dlqi_absolute <- function(df) {
  b <- baseline_lookup(df, "abs_dlqi_mean", "abs_dlqi_sd")
  d <- derive_change(b$mean, df$abs_dlqi_mean, df$abs_dlqi_change_mean)
  df$baseline        <- fmt_mean_sd(b$mean, b$sd)
  df$on_tx           <- fmt_mean_sd(df$abs_dlqi_mean, df$abs_dlqi_sd)
  df$abs_dlqi_change <- fmt_mean_sd_marked(d$change, df$abs_dlqi_change_sd,
                                           d$change_derived)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  df <- df[is.na(df$timepoint) | df$timepoint > 0, ]
  has_any <- nzchar(df$baseline) | nzchar(df$on_tx) | nzchar(df$abs_dlqi_change)
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", "baseline", "on_tx", "abs_dlqi_change")]
}

# Endpoint catalogue. Each tab has one or more endpoint groups. A group is
# defined by the source v_* table, the format fn that shapes it for display,
# and the column headers shown to the user.
endpoint_groups <- list(
  pasi = list(
    label = "PASI",
    groups = list(
      response = list(
        label    = "PASI 50 / 75 / 90 / 100",
        table    = "v_pasi",
        fmt      = format_pasi_response,
        colnames = c("Trial", "Drug", "N",
                     "PASI 50", "PASI 75", "PASI 90", "PASI 100")
      ),
      absolute = list(
        label    = "Absolute PASI",
        table    = "v_pasi_abs",
        fmt      = format_pasi_absolute,
        colnames = c("Trial", "Drug", "N",
                     "Baseline", "Follow-up", "Δ from baseline"),
        note     = paste("Change values marked with * are derived as",
                         "follow-up - baseline when the study reported",
                         "both timepoints but not the change directly;",
                         "SDs are not shown for derived values.")
      )
    )
  ),
  dlqi = list(
    label = "DLQI",
    groups = list(
      zero = list(
        label    = "DLQI 0/1, DLQI 0",
        table    = "v_dlqi",
        fmt      = format_dlqi_zero,
        colnames = c("Trial", "Drug", "N", "DLQI 0/1", "DLQI 0")
      ),
      threshold = list(
        label    = "DLQI ≤ 5",
        table    = "v_dlqi",
        fmt      = format_dlqi_threshold,
        colnames = c("Trial", "Drug", "N", "DLQI ≤ 5")
      ),
      change = list(
        label    = "5+ / 4+ point decrease",
        table    = "v_dlqi",
        fmt      = format_dlqi_change,
        colnames = c("Trial", "Drug", "N",
                     "5+ pt decrease", "4+ pt decrease")
      ),
      absolute = list(
        label    = "Absolute DLQI",
        table    = "v_dlqi",
        fmt      = format_dlqi_absolute,
        colnames = c("Trial", "Drug", "N",
                     "Baseline", "Follow-up", "Δ from baseline"),
        note     = paste("Change values marked with * are derived as",
                         "follow-up - baseline when the study reported",
                         "both timepoints but not the change directly;",
                         "SDs are not shown for derived values.")
      )
    )
  ),
  safety = list(
    label = "Safety",
    groups = list(
      sae = list(
        label    = "Any SAE",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df, c("sae")),
        colnames = c("Trial", "Drug", "N", "Any SAE")
      ),
      disc = list(
        label    = "Discontinuation",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df,
                     c("disc_any", "disc_ae")),
        colnames = c("Trial", "Drug", "N",
                     "Disc. (any)", "Disc. (AE)")
      ),
      serious_infection = list(
        label    = "Serious infections",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df, c("serious_infection")),
        colnames = c("Trial", "Drug", "N", "Serious infection")
      ),
      injection_site_rxn = list(
        label    = "Injection-site reactions",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df, c("injection_site_rxn")),
        colnames = c("Trial", "Drug", "N", "Injection site rxn")
      ),
      malignancy = list(
        label    = "Malignancies",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df,
                     c("malignancy", "nmsc", "malignancy_non_nmsc")),
        colnames = c("Trial", "Drug", "N",
                     "Malignancy", "NMSC", "Malignancy (non-NMSC)")
      )
    )
  )
)

# Per-group "does this (trial, drug) arm contribute a row to this endpoint?"
# predicate. Mirrors the row-drop logic inside each fmt() in endpoint_groups
# (e.g. `has_any` across the binary cols) but also computes a per-arm
# patient count (max n across the arm's surviving rows in this endpoint,
# i.e. the ITT denominator for the endpoint). Returned columns drive both
# the per-endpoint network subset and node sizing.
empty_survivors <- function() {
  data.frame(trial = character(), ref_id = integer(), arm_no = integer(),
             drug = character(), n_arm = integer(),
             stringsAsFactors = FALSE)
}

summarise_arms <- function(df) {
  if (!nrow(df)) return(empty_survivors())
  agg <- aggregate(df$n,
                   by = list(trial = df$trial, ref_id = df$ref_id,
                             arm_no = df$arm_no, drug = df$drug),
                   FUN = function(x) {
                     v <- suppressWarnings(max(x, na.rm = TRUE))
                     if (is.infinite(v)) NA_integer_ else as.integer(v)
                   })
  names(agg)[ncol(agg)] <- "n_arm"
  agg
}

survivors_binary <- function(df, cols) {
  if (!nrow(df)) return(empty_survivors())
  ok <- Reduce(`|`, lapply(cols, function(c)
    !is.na(df[[c]]) & !is.na(df$n) & df$n > 0))
  summarise_arms(df[ok, , drop = FALSE])
}

# Absolute-value endpoints: a (trial, drug) arm contributes iff some
# non-baseline timepoint has a mean or change-from-baseline value. Matches
# format_pasi_absolute / format_dlqi_absolute, which drop t==0 rows and then
# require any of {baseline, on_tx, change} to be non-empty. Baseline alone
# without a follow-up value would not survive their `has_any` either.
survivors_absolute <- function(df, mean_col, change_col) {
  if (!nrow(df)) return(empty_survivors())
  tp_ok   <- is.na(df$timepoint) | df$timepoint > 0
  any_val <- !is.na(df[[mean_col]]) | !is.na(df[[change_col]])
  summarise_arms(df[tp_ok & any_val, , drop = FALSE])
}

group_survivors <- list(
  pasi = list(
    response = function(df) survivors_binary(df,
                              c("pasi50", "pasi75", "pasi90", "pasi100")),
    absolute = function(df) survivors_absolute(df,
                              "abs_pasi_mean", "abs_pasi_change_mean")
  ),
  dlqi = list(
    zero      = function(df) survivors_binary(df,  c("dlqi_0_1", "dlqi_0")),
    threshold = function(df) survivors_binary(df,  c("dlqi_le5")),
    change    = function(df) survivors_binary(df,
                               c("dlqi_5pt_dec", "dlqi_4pt_dec")),
    absolute  = function(df) survivors_absolute(df,
                               "abs_dlqi_mean", "abs_dlqi_change_mean")
  ),
  safety = list(
    sae                = function(df) survivors_binary(df, c("sae")),
    disc               = function(df) survivors_binary(df,
                                        c("disc_any", "disc_ae")),
    serious_infection  = function(df) survivors_binary(df,
                                        c("serious_infection")),
    injection_site_rxn = function(df) survivors_binary(df,
                                        c("injection_site_rxn")),
    malignancy         = function(df) survivors_binary(df,
                                        c("malignancy", "nmsc",
                                          "malignancy_non_nmsc"))
  )
)

# Precompute the (trial, drug) set and the network for each (tab, group).
# Switching endpoints is then a free dictionary lookup at runtime.
group_td <- list()
for (.tid in names(endpoint_groups)) {
  group_td[[.tid]] <- list()
  for (.gid in names(endpoint_groups[[.tid]]$groups)) {
    .grp <- endpoint_groups[[.tid]]$groups[[.gid]]
    .raw <- read_db(sprintf("SELECT * FROM %s", .grp$table))
    group_td[[.tid]][[.gid]] <- group_survivors[[.tid]][[.gid]](.raw)
  }
}

# Find the largest per-drug patient total across any single endpoint group.
# This anchors node sizing so the same patient count produces the same
# visual size in every subset (and the biggest drug in the busiest endpoint
# hits NODE_SIZE_MAX exactly).
.global_max_patients <- max(0L, unlist(lapply(group_td, function(tab_grps) {
  unlist(lapply(tab_grps, function(td) {
    if (!nrow(td)) return(0L)
    td <- unique(td[, c("trial", "ref_id", "arm_no", "drug", "n_arm")])
    if (!nrow(td)) return(0L)
    vapply(unique(td$drug), function(d)
      sum(td$n_arm[td$drug == d], na.rm = TRUE), numeric(1))
  }))
})), na.rm = TRUE)

# Master network = union of every endpoint group's (trial, drug) set. Its
# layout is computed once with layout_with_kk and reused for every subset,
# so a drug's position never moves when the user switches endpoints.
master_td <- unique(do.call(rbind, lapply(group_td, function(tab_grps)
  do.call(rbind, tab_grps))))
master_network <- build_network_data(master_td, ref_max_n = .global_max_patients)

.master_layout <- visNetwork(master_network$nodes, master_network$edges) |>
  visIgraphLayout(layout = "layout_with_kk", randomSeed = 42, physics = FALSE)
.layout_nodes <- .master_layout$x$nodes
.mi <- match(master_network$nodes$id, .layout_nodes$id)
# visIgraphLayout normalises coords to ~[-1, 1] and relies on a JS-side
# square-fit multiplier to spread them across the canvas. We bypass that
# wrapper (so subset switches don't relayout), so scale to pixel space
# ourselves — ~650 gives the larger patient-count-sized nodes room to
# breathe without overlapping labels.
master_network$nodes$x <- .layout_nodes$x[.mi] * 650
master_network$nodes$y <- .layout_nodes$y[.mi] * 650

attach_master_coords <- function(nw) {
  m <- match(nw$nodes$id, master_network$nodes$id)
  nw$nodes$x <- master_network$nodes$x[m]
  nw$nodes$y <- master_network$nodes$y[m]
  nw
}

group_networks <- list()
for (.tid in names(endpoint_groups)) {
  group_networks[[.tid]] <- list()
  for (.gid in names(endpoint_groups[[.tid]]$groups)) {
    group_networks[[.tid]][[.gid]] <- attach_master_coords(
      build_network_data(group_td[[.tid]][[.gid]],
                         ref_max_n = .global_max_patients)
    )
  }
}

# Preserved for downstream lookups (filter-state click handlers, etc.).
nodes_df <- master_network$nodes
edges_df <- master_network$edges

# ---------------------------------------------------------------------------
# Meta-analysis catalogue. Maps each (tab_id, group_id) in `endpoint_groups`
# to the precomputed MA outcomes drawn in the modal. `scale` controls the
# x-axis treatment ("rr" = log scale, ref line at 1; "md" = linear, ref 0;
# "prop" = 0..1, no ref line). `label` is the plot title.
# ---------------------------------------------------------------------------
ma_catalog <- list(
  pasi = list(
    response = list(
      list(code = "pasi50",  label = "PASI 50", scale = "rr"),
      list(code = "pasi75",  label = "PASI 75", scale = "rr"),
      list(code = "pasi90",  label = "PASI 90", scale = "rr"),
      list(code = "pasi100", label = "PASI 100", scale = "rr")
    ),
    absolute = list(
      list(code = "abs_pasi_change", label = "Δ from baseline (absolute PASI)",
           scale = "md", endpoint_group = "pasi_abs")
    )
  ),
  dlqi = list(
    zero = list(
      list(code = "dlqi_0_1", label = "DLQI 0/1", scale = "rr"),
      list(code = "dlqi_0",   label = "DLQI 0",   scale = "rr")
    ),
    threshold = list(
      list(code = "dlqi_le5", label = "DLQI ≤ 5", scale = "rr")
    ),
    change = list(
      list(code = "dlqi_5pt_dec", label = "5+ point decrease", scale = "rr"),
      list(code = "dlqi_4pt_dec", label = "4+ point decrease", scale = "rr")
    ),
    absolute = list(
      list(code = "abs_dlqi_change", label = "Δ from baseline (absolute DLQI)",
           scale = "md", endpoint_group = "dlqi_abs")
    )
  ),
  safety = list(
    sae                = list(list(code = "sae", label = "Any SAE", scale = "rr")),
    disc               = list(
      list(code = "disc_any", label = "Discontinuation (any)", scale = "rr"),
      list(code = "disc_ae",  label = "Discontinuation (AE)", scale = "rr")
    ),
    serious_infection  = list(
      list(code = "serious_infection", label = "Serious infection", scale = "rr")
    ),
    injection_site_rxn = list(
      list(code = "injection_site_rxn", label = "Injection-site reaction",
           scale = "rr")
    ),
    malignancy = list(
      list(code = "malignancy",          label = "Malignancy", scale = "rr"),
      list(code = "nmsc",                label = "NMSC", scale = "rr"),
      list(code = "malignancy_non_nmsc", label = "Malignancy (non-NMSC)",
           scale = "rr")
    )
  )
)

# Endpoint-group key used inside the MA tables. Defaults to the tab id,
# overridden by `endpoint_group` on outcomes that live in a different bucket
# (the absolute-value endpoints).
ma_endpoint_group <- function(tab_id, outcome) {
  outcome$endpoint_group %||% tab_id
}

# Single source of truth for whether the MA build step has been run.
ma_tables_present <- function() {
  con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  tabs <- dbListTables(con)
  all(c("ma_pairwise", "ma_pairwise_trials",
        "ma_proportion", "ma_proportion_trials") %in% tabs)
}
HAS_MA <- ma_tables_present()
MA_BUILT_AT <- if (HAS_MA) {
  tryCatch(read_db("SELECT MAX(built_at) AS b FROM ma_pairwise")$b,
           error = function(e) NA_character_)
} else NA_character_

# ---------------------------------------------------------------------------
# Forest-plot helpers. One renderer (`forest_ggiraph`) is reused for all
# three modal kinds; the only thing that varies is what each row represents
# (per-trial squares for "pairwise" and "proportion", per-drug squares for
# "drug_vs_placebo"). Tooltips carry trial/drug name, n / N, effect, CI
# and FE+RE weights.
# ---------------------------------------------------------------------------
fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f (%.*f to %.*f)",
          digits, est, digits, lo, digits, hi)
}

fmt_pct <- function(x, digits = 0) {
  sprintf("%.*f%%", digits, 100 * x)
}

# Build a tooltip string (HTML) for a forest-plot row. `extra` is an
# optional named character vector of additional label/value pairs.
ma_tooltip <- function(label, est, lo, hi, weight_fe, weight_re, extra = NULL,
                       digits = 2) {
  parts <- c(
    sprintf("<b>%s</b>", htmltools::htmlEscape(label)),
    sprintf("Effect: %s", fmt_ci(est, lo, hi, digits)),
    if (!is.na(weight_fe))
      sprintf("Weight (FE): %s", fmt_pct(weight_fe, 1)),
    if (!is.na(weight_re))
      sprintf("Weight (RE): %s", fmt_pct(weight_re, 1))
  )
  if (length(extra))
    parts <- c(parts, sprintf("%s: %s", names(extra), as.character(extra)))
  paste(parts, collapse = "<br/>")
}

# Render one forest plot. `rows` is a data.frame with columns:
#   label, est, lo, hi, weight_fe, weight_re, tooltip
# Optional `pooled` rows are drawn as diamonds at the bottom; each pooled
# row has label, est, lo, hi, kind ("FE" or "RE"). `scale` is one of
# "rr", "md", "prop". `emphasis` is "FE" or "RE" (the highlighted diamond).
forest_ggiraph <- function(rows, pooled, scale = "rr", emphasis = "RE",
                            title = NULL, output_id) {
  if (!nrow(rows)) {
    return(girafe(ggobj = ggplot() +
                    annotate("text", x = 0, y = 0, label = "No data") +
                    theme_void(),
                  width_svg = 7, height_svg = 1.5))
  }

  ref_line <- switch(scale, rr = 1, md = 0, prop = NA_real_)
  x_log    <- identical(scale, "rr")

  # Order rows from top (first) to bottom (last); ggplot's y is numeric
  # going up, so reverse.
  rows$y <- rev(seq_len(nrow(rows)))
  rows$size_w <- if (all(is.na(rows$weight_fe))) 0.5
                 else pmax(0.15, sqrt(rows$weight_fe / max(rows$weight_fe, na.rm = TRUE)))

  # Diamonds for pooled rows, vertically below all trial rows.
  pooled_geom <- NULL
  if (!is.null(pooled) && nrow(pooled)) {
    pooled$y <- -seq_len(nrow(pooled)) * 0.8
    diamonds <- do.call(rbind, lapply(seq_len(nrow(pooled)), function(i) {
      p <- pooled[i, ]
      hh <- 0.28
      data.frame(
        kind = p$kind,
        x = c(p$lo, p$est, p$hi, p$est),
        y = p$y + c(0, hh, 0, -hh)
      )
    }))
    diamonds$y_id <- factor(rep(pooled$kind, each = 4))
    pooled_geom <- list(
      geom_polygon(data = diamonds,
                   aes(x = x, y = y, group = y_id, fill = kind, colour = kind),
                   alpha = 0.55),
      geom_text(data = pooled,
                aes(x = max(rows$hi, hi, na.rm = TRUE) * (if (x_log) 1.6 else 1),
                    y = y,
                    label = sprintf("%s: %s", kind,
                                    fmt_ci(est, lo, hi,
                                           if (scale == "prop") 3 else 2))),
                hjust = 0, size = 3, colour = "#1a1f2c")
    )
  }

  fill_vals <- c(FE = "#7f8fa6", RE = "#1F4E8C")
  if (emphasis == "FE") fill_vals <- c(FE = "#1F4E8C", RE = "#c9d3df")
  else                  fill_vals <- c(FE = "#c9d3df", RE = "#1F4E8C")

  p <- ggplot(rows, aes(x = est, y = y)) +
    geom_vline(xintercept = ref_line, linetype = "dashed",
               colour = "#aab1bd", na.rm = TRUE) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0,
                  orientation = "y", colour = "#5a6478", na.rm = TRUE) +
    geom_point_interactive(aes(size = size_w, tooltip = tooltip,
                                data_id = as.character(y)),
                            shape = 15, colour = "#1a1f2c", na.rm = TRUE) +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_y_continuous(breaks = rows$y, labels = rows$label,
                       expand = expansion(add = c(2.5, 0.6))) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(colour = "#1a1f2c"),
          legend.position = "top",
          plot.title = element_text(size = 12, face = "bold",
                                    colour = "#1a1f2c"))
  if (x_log) p <- p + scale_x_log10()
  if (!is.null(pooled_geom)) {
    p <- p + pooled_geom +
      scale_fill_manual(values = fill_vals, name = NULL) +
      scale_colour_manual(values = fill_vals, guide = "none")
  }

  height_in <- 0.5 + 0.25 * nrow(rows) +
               (if (!is.null(pooled)) 0.4 * nrow(pooled) else 0)
  girafe(ggobj = p, width_svg = 9, height_svg = max(2.2, height_in),
         options = list(
           opts_tooltip(css = "background:#1a1f2c;color:#fff;padding:6px 9px;
                              border-radius:5px;font-size:12px;
                              box-shadow:0 2px 8px rgba(0,0,0,0.25);"),
           opts_hover(css = "stroke:#c0392b;stroke-width:2px;"),
           opts_sizing(rescale = TRUE, width = 1)
         ))
}

# --- Data fetchers used by the modal ---------------------------------------

# Pull a single pairwise MA row + its trial rows. Returns NULL if absent.
fetch_pairwise <- function(drug_a, drug_b, endpoint_group, outcome_code) {
  r <- read_db(
    "SELECT * FROM ma_pairwise
     WHERE drug_a = ? AND drug_b = ?
       AND endpoint_group = ? AND outcome_code = ?",
    params = list(drug_a, drug_b, endpoint_group, outcome_code))
  if (!nrow(r)) return(NULL)
  trials <- read_db(
    "SELECT * FROM ma_pairwise_trials WHERE comparison_id = ?",
    params = list(r$comparison_id[1]))
  list(summary = r[1, , drop = FALSE], trials = trials)
}

# All drug-vs-placebo MA rows for one outcome (no per-trial detail).
fetch_drugs_vs_placebo <- function(endpoint_group, outcome_code) {
  read_db(
    "SELECT * FROM ma_pairwise
     WHERE drug_b = 'Placebo' AND endpoint_group = ? AND outcome_code = ?
     ORDER BY drug_a",
    params = list(endpoint_group, outcome_code))
}

# Single-arm proportion + per-trial rows for one (drug, outcome).
fetch_proportion <- function(drug, endpoint_group, outcome_code) {
  r <- read_db(
    "SELECT * FROM ma_proportion
     WHERE drug = ? AND endpoint_group = ? AND outcome_code = ?",
    params = list(drug, endpoint_group, outcome_code))
  if (!nrow(r)) return(NULL)
  trials <- read_db(
    "SELECT * FROM ma_proportion_trials WHERE proportion_id = ?",
    params = list(r$proportion_id[1]))
  list(summary = r[1, , drop = FALSE], trials = trials)
}

# Translate one MA "outcome spec" + the current filter into a rows/pooled
# bundle for `forest_ggiraph`. Returns NULL when no data is available.
build_forest_inputs <- function(state, tab_id, outcome) {
  endpoint_group <- ma_endpoint_group(tab_id, outcome)
  digits <- if (identical(outcome$scale, "md")) 2 else 2
  if (is.null(state)) {
    # No filter -> one row per drug vs placebo.
    df <- fetch_drugs_vs_placebo(endpoint_group, outcome$code)
    if (!nrow(df)) return(NULL)
    rows <- data.frame(
      label = df$drug_a,
      est   = df$te_re,
      lo    = df$lo_re,
      hi    = df$hi_re,
      weight_fe = NA_real_,
      weight_re = NA_real_,
      stringsAsFactors = FALSE
    )
    if (identical(outcome$scale, "rr")) {
      rows$est <- exp(rows$est); rows$lo <- exp(rows$lo); rows$hi <- exp(rows$hi)
    }
    rows$tooltip <- vapply(seq_len(nrow(rows)), function(i)
      ma_tooltip(rows$label[i], rows$est[i], rows$lo[i], rows$hi[i],
                 NA_real_, NA_real_,
                 extra = c("n trials" = df$n_studies[i],
                           "I² " = sprintf("%.1f%%", 100 * df$i2[i]))),
      character(1))
    return(list(rows = rows, pooled = NULL))
  }
  if (identical(state$kind, "edge")) {
    # Edge filter -> per-trial pairwise forest, with FE + RE diamonds.
    # Order (a, b) so the stored direction is hit: drug_a vs drug_b. We
    # check both orientations to be tolerant.
    for (orient in list(c(state$from, state$to), c(state$to, state$from))) {
      res <- fetch_pairwise(orient[1], orient[2], endpoint_group, outcome$code)
      if (!is.null(res)) break
    }
    if (is.null(res)) return(NULL)
    s  <- res$summary; t <- res$trials
    rows <- data.frame(
      label = t$trial,
      est   = t$te, lo = t$lo, hi = t$hi,
      weight_fe = t$weight_fe, weight_re = t$weight_re,
      stringsAsFactors = FALSE
    )
    pooled <- data.frame(
      kind = c("FE", "RE"),
      est  = c(s$te_fe, s$te_re),
      lo   = c(s$lo_fe, s$lo_re),
      hi   = c(s$hi_fe, s$hi_re),
      stringsAsFactors = FALSE
    )
    if (identical(outcome$scale, "rr")) {
      rows$est <- exp(rows$est); rows$lo <- exp(rows$lo); rows$hi <- exp(rows$hi)
      pooled$est <- exp(pooled$est); pooled$lo <- exp(pooled$lo); pooled$hi <- exp(pooled$hi)
    }
    rows$tooltip <- vapply(seq_len(nrow(rows)), function(i) {
      extra <- if (identical(outcome$scale, "rr"))
        c("Arm A" = sprintf("%d / %d", t$event_a[i], t$n_a[i]),
          "Arm B" = sprintf("%d / %d", t$event_b[i], t$n_b[i]))
      else
        c("A: mean (SD), n" = sprintf("%.2f (%.2f), %d", t$mean_a[i], t$sd_a[i], t$n_a[i]),
          "B: mean (SD), n" = sprintf("%.2f (%.2f), %d", t$mean_b[i], t$sd_b[i], t$n_b[i]))
      ma_tooltip(rows$label[i], rows$est[i], rows$lo[i], rows$hi[i],
                 rows$weight_fe[i], rows$weight_re[i],
                 extra = extra)
    }, character(1))
    return(list(rows = rows, pooled = pooled, comparison = sprintf("%s vs %s",
                                                                    s$drug_a, s$drug_b)))
  }
  if (identical(state$kind, "node")) {
    # Node filter -> single-arm proportion (binary outcomes only).
    if (!identical(outcome$scale, "rr")) return(NULL)
    res <- fetch_proportion(state$drug, endpoint_group, outcome$code)
    if (is.null(res)) return(NULL)
    s <- res$summary; t <- res$trials
    rows <- data.frame(
      label = t$trial,
      est   = t$p, lo = t$lo, hi = t$hi,
      weight_fe = t$weight_fe, weight_re = t$weight_re,
      stringsAsFactors = FALSE
    )
    rows$tooltip <- vapply(seq_len(nrow(rows)), function(i)
      ma_tooltip(rows$label[i], rows$est[i], rows$lo[i], rows$hi[i],
                 rows$weight_fe[i], rows$weight_re[i],
                 extra = c("Responders" = sprintf("%d / %d",
                                                  t$k[i], t$n[i])),
                 digits = 3),
      character(1))
    pooled <- data.frame(
      kind = c("FE", "RE"),
      est  = c(s$te_fe, s$te_re),
      lo   = c(s$lo_fe, s$lo_re),
      hi   = c(s$hi_fe, s$hi_re),
      stringsAsFactors = FALSE
    )
    # Proportions are already back-transformed in storage; force linear scale.
    return(list(rows = rows, pooled = pooled, drug = state$drug,
                effective_scale = "prop"))
  }
  NULL
}

ui <- fluidPage(
  tags$head(tags$script(HTML("
    // Disable specific <option> values inside a Shiny select input. Used to
    // grey out endpoint groups that have zero rows under the current filter.
    Shiny.addCustomMessageHandler('set_disabled_options', function(msg) {
      var sel = document.getElementById(msg.input_id);
      if (!sel) return;
      var disabled = Array.isArray(msg.disabled) ? msg.disabled : [];
      Array.from(sel.options).forEach(function(opt) {
        opt.disabled = disabled.indexOf(opt.value) !== -1;
      });
    });

    // Trial-name popover. fmt_trial() emits <a class='trial-pop'
    // data-citations='...'> elements; clicking one opens a floating panel
    // anchored beneath the link. Event delegation on document means DT
    // redraws (filter/paginate) don't break the trigger.
    (function() {
      var pop = null;
      function ensure() {
        if (pop) return pop;
        pop = document.createElement('div');
        pop.id = 'trial-popover';
        pop.style.display = 'none';
        document.body.appendChild(pop);
        return pop;
      }
      function hide() { if (pop) pop.style.display = 'none'; }
      function show(trigger) {
        var p = ensure();
        var id = trigger.getAttribute('data-ref-id');
        p.innerHTML = (window.studyCitations && window.studyCitations[id]) || '';
        p.style.display = 'block';
        // Measure after layout so we can flip when there's no room below.
        var r  = trigger.getBoundingClientRect();
        var ph = p.offsetHeight;
        var pw = p.offsetWidth;
        var spaceBelow = window.innerHeight - r.bottom;
        var top = (spaceBelow < ph + 12 && r.top > spaceBelow)
          ? window.scrollY + r.top - ph - 6
          : window.scrollY + r.bottom + 6;
        var left = window.scrollX + r.left;
        var maxLeft = window.scrollX + window.innerWidth - pw - 12;
        if (left > maxLeft) left = Math.max(8, maxLeft);
        p.style.top  = top  + 'px';
        p.style.left = left + 'px';
      }
      document.addEventListener('click', function(ev) {
        var t = ev.target.closest && ev.target.closest('.trial-pop');
        if (t) { ev.preventDefault(); show(t); return; }
        if (pop && pop.contains(ev.target)) return;  // clicks inside popover
        hide();
      });
      document.addEventListener('keydown', function(ev) {
        if (ev.key === 'Escape') hide();
      });
      window.addEventListener('scroll', hide, true);
    })();
  "))),
  tags$head(tags$script(HTML(paste0(
    "window.studyCitations = ", cites_json, ";"
  )))),
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com",
              crossorigin = NA),
    tags$link(rel = "stylesheet",
              href = paste0("https://fonts.googleapis.com/css2?",
                            "family=Inter:wght@400;500;600;700&display=swap"))
  ),
  tags$head(tags$style(HTML("
    /* --- base typography ------------------------------------------------ */
    body, .container-fluid, .container, button, input, select, textarea {
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI',
                   Roboto, Helvetica, Arial, sans-serif;
      color: #1a1f2c;
      -webkit-font-smoothing: antialiased;
    }
    body { background: #ffffff; color: #1a1f2c; }
    h1, h2, h3, h4, h5 { font-weight: 600; letter-spacing: -0.01em; }
    /* App header: title left, action buttons right, sits as a unified bar
       across the top of the page with a subtle separator below. Negative
       horizontal margins bleed it past the container's 15px gutters so
       the bottom border runs edge-to-edge. */
    .title-bar {
      background: #ffffff;
      border-bottom: 1px solid #e3e6ea;
      padding: 10px 24px;
      margin: 0 -15px 12px -15px;
      flex: 0 0 auto;
    }
    /* Trim the default Bootstrap titlePanel — it's huge and bold. */
    .title-bar h2 { font-size: 22px; font-weight: 600; margin: 0;
                    letter-spacing: -0.015em; color: #0f1726;
                    line-height: 1.2; }
    .title-bar .title-actions .btn { padding: 5px 12px; }

    /* --- buttons -------------------------------------------------------- */
    .btn-default {
      background: #ffffff; border: 1px solid #d4d8df; color: #2a3142;
      border-radius: 5px; font-weight: 500; transition: all 0.15s ease;
      box-shadow: none;
    }
    .btn-default:hover, .btn-default:focus {
      background: #f5f7fa; border-color: #b8bec8; color: #0f1726;
    }
    .btn-default:active { background: #eceff3; }

    /* --- tabs: flat underline style ------------------------------------- */
    .nav-tabs {
      border-bottom: 1px solid #e3e6ea; margin-bottom: 0; padding-left: 2px;
    }
    .nav-tabs > li { margin-bottom: -1px; }
    .nav-tabs > li > a {
      border: none; border-bottom: 2px solid transparent;
      background: transparent; color: #5a6478; font-weight: 500;
      padding: 10px 16px; margin-right: 4px; border-radius: 0;
      transition: color 0.15s ease, border-color 0.15s ease;
    }
    .nav-tabs > li > a:hover {
      background: transparent; border-bottom-color: #c8ced6; color: #1F4E8C;
    }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:hover,
    .nav-tabs > li.active > a:focus {
      background: transparent; border: none;
      border-bottom: 2px solid #1F4E8C; color: #1F4E8C; font-weight: 600;
    }

    .title-actions { display: flex; align-items: center; gap: 8px; }
    .title-actions .btn { padding: 4px 12px; font-size: 13px; }
    .view-summary-filter { font-weight: 600; color: #1F4E8C; }
    .view-summary-sep    { color: #5a6478; }
    #nma { background: #fafbfc; border: 1px solid #e3e6ea; border-radius: 6px; }
    .endpoint-picker { margin: 0; padding: 6px 0; background: #ffffff;
                       height: 46px; box-sizing: border-box; }
    .endpoint-picker .form-group { margin-bottom: 0; }
    .endpoint-picker select { width: 100%; height: 34px; padding: 4px 8px;
                              border: 1px solid #ccc; border-radius: 4px;
                              background: #ffffff; font-size: 14px; }
    .endpoint-picker select option:disabled { color: #b0b0b0; }
    /* Use flexbox so left/right columns share the row's full height; sticky
       then has room to stick as the user scrolls the table. */
    /* Lock the page so only the table column scrolls; the diagram never moves. */
    html, body { height: 100%; overflow: hidden; }
    .container-fluid, .container { height: 100%; display: flex; flex-direction: column; }
    .row.split { display: flex; align-items: stretch; flex: 1 1 auto;
                 min-height: 0; }
    .row.split > [class*='col-'] { float: none; }
    /* Only the table body scrolls; tabs stay pinned at the top of the column. */
    .row.split .col-table { display: flex; flex-direction: column;
                            max-height: 100%; min-height: 0; }
    .row.split .col-table .tabbable { display: flex; flex-direction: column;
                                       flex: 1 1 auto; min-height: 0; }
    .row.split .col-table .tab-content { flex: 1 1 auto; overflow: auto;
                                          min-height: 0; }
    /* Pin the endpoint dropdown at the top of the scrolling tab-content,
       and pin the summary line + table header directly underneath it. */
    .row.split .col-table .tab-content .endpoint-picker {
      position: sticky; top: 0; z-index: 3;
    }
    /* Let the sticky .view-summary's containing block be the scrolling
       .tab-content rather than this Shiny wrapper (which has no extra
       height of its own and would pin sticky in place). */
    .view-summary-wrap { display: contents; }
    .row.split .col-table .tab-content .view-summary {
      position: sticky; top: 46px; z-index: 3;
      background: #ffffff; padding: 10px 4px 10px 4px;
      font-size: 16px; font-weight: 500; color: #1a1f2c;
      border-bottom: 1px solid #e3e6ea;
    }
    .row.split .col-table .tab-content table.dataTable thead th {
      position: sticky; top: 90px; background: #ffffff; z-index: 2;
      box-shadow: inset 0 -1px 0 #ddd;
    }
    .app-footer { margin-top: 12px; font-size: 12px; color: #888; }
    .app-footer a { color: #1F4E8C; text-decoration: none; }
    .app-footer a:hover { text-decoration: underline; }
    /* Trial-name popover. Anchored absolutely; positioning is set inline by
       the click handler in <script> above. */
    a.trial-pop { color: #1F4E8C; cursor: pointer; }
    a.trial-pop:hover { text-decoration: underline; }
    #trial-popover {
      position: absolute; z-index: 1000; max-width: 520px;
      background: #ffffff; border: 1px solid #c8ced6; border-radius: 6px;
      box-shadow: 0 4px 16px rgba(0,0,0,0.12);
      padding: 10px 14px; font-size: 13px; line-height: 1.45; color: #222;
    }
    #trial-popover .trial-cite { margin-bottom: 8px; }
    #trial-popover .trial-cite:last-child { margin-bottom: 0; }
    #trial-popover a { color: #1F4E8C; word-break: break-all; }
    /* Per-trial banding applied by drawCallback above. */
    table.dataTable tr.trial-band > td { background-color: #f3f5f8; }
    .title-bar { display: flex; align-items: center; justify-content: space-between;
                 gap: 12px; }
    .title-bar .about-btn { margin-right: 4px; }
    /* Roomy About modal. */
    .about-modal .modal-dialog { width: 80%; max-width: 900px; }
    .about-modal .modal-body { font-size: 14px; line-height: 1.55;
                               max-height: 70vh; overflow-y: auto; }
    .about-modal .modal-body h4 { margin-top: 18px; color: #1F4E8C; }
    .about-modal .modal-body h4:first-child { margin-top: 0; }
    /* Meta-analysis modal: a touch wider than About for forest plots. */
    .ma-modal .modal-dialog { width: 92%; max-width: 1100px; }
    .ma-modal .modal-body { max-height: 75vh; overflow-y: auto;
                            padding: 14px 18px; }
    .ma-modal h4 { margin-top: 4px; color: #1F4E8C; font-size: 15px; }
    .ma-modal .ma-summary { font-size: 13px; color: #5a6478;
                            margin-bottom: 14px; }
    .ma-modal .ma-plot-block { margin-bottom: 22px; }
    .ma-modal .ma-plot-title { font-weight: 600; color: #1a1f2c;
                               font-size: 14px; margin: 8px 0 2px 0; }
    .ma-modal .ma-plot-stats { font-size: 12px; color: #5a6478;
                               margin-bottom: 6px; }
    .ma-modal .ma-empty { color: #5a6478; font-style: italic;
                          padding: 12px 0; }
    .ma-modal .ma-footer-meta { font-size: 12px; color: #888;
                                margin-right: auto; }
    .ma-modal .ma-toggle { display: inline-flex; align-items: center;
                            gap: 6px; margin-right: 12px; font-size: 13px; }
    .ma-modal .ma-toggle .radio-inline { margin-left: 6px; }
  "))),
  div(class = "title-bar",
      titlePanel("Psoriasis RCT Explorer"),
      div(class = "title-actions",
          downloadButton("download_db", "Download SQLite",
                         class = "btn btn-default"),
          actionButton("show_ma", "Meta-analyse",
                       icon = icon("chart-column"),
                       class = "btn btn-default ma-btn"),
          actionButton("show_about", "About", icon = icon("circle-info"),
                       class = "btn btn-default about-btn"))
  ),
  fluidRow(class = "split",
    column(6,
      visNetworkOutput("nma", height = "640px"),
      helpText("Click a node to filter to one drug; click an edge to show only",
               "trials comparing that pair. The network reflects the currently",
               "selected endpoint — drugs and comparisons without data for that",
               "endpoint are hidden. Node area is proportional to the number",
               "of randomised patients contributing to the endpoint; edge",
               "width is the number of trials with that head-to-head. Click",
               "empty space to clear the filter."),
      tags$footer(class = "app-footer",
        HTML("&copy; 2026 Thomas Clemmet"))
    ),
    column(6, class = "col-table",
      do.call(tabsetPanel, c(
        list(id = "view", type = "tabs"),
        lapply(names(endpoint_groups), function(tab_id) {
          tab <- endpoint_groups[[tab_id]]
          group_choices <- setNames(names(tab$groups),
                                    vapply(tab$groups, `[[`, "", "label"))
          tabPanel(
            tab$label,
            value = tab_id,
            div(class = "endpoint-picker",
                selectInput(paste0("group_", tab_id),
                            label = NULL,
                            choices  = group_choices,
                            selected = group_choices[[1]],
                            selectize = FALSE,
                            width    = "100%")),
            uiOutput(paste0("summary_", tab_id),
                     class = "view-summary-wrap"),
            DTOutput(paste0("tbl_", tab_id))
          )
        })
      ))
    )
  )
)

server <- function(input, output, session) {

  filter_state <- reactiveVal(NULL)

  # Plain-text rendering of the current filter for the per-tab summary line.
  filter_text <- function(s) {
    if (is.null(s))            "All drugs"
    else if (s$kind == "node") s$drug
    else if (s$kind == "edge") sprintf("%s ↔ %s", s$from, s$to)
    else                       "All drugs"
  }

  # For each (tab, group), does the current filter yield any rows after
  # formatting? Cache table queries within one pass since v_safety is reused
  # across five groups.
  availability <- reactive({
    state <- filter_state()
    cache <- list()
    get_tbl <- function(tbl) {
      if (is.null(cache[[tbl]])) cache[[tbl]] <<- query_view(tbl, state)
      cache[[tbl]]
    }
    lapply(endpoint_groups, function(tab) {
      vapply(tab$groups, function(grp) {
        n <- tryCatch(nrow(grp$fmt(get_tbl(grp$table))),
                      error = function(e) 0L)
        isTRUE(n > 0)
      }, logical(1))
    })
  })

  # Push disabled-option lists into each tab's native <select>. We don't
  # auto-switch the selection: that would re-render the same DT widget with
  # a different column count and DataTables.js fires a "column not found"
  # warning. Keeping the user's pick lets the empty-table path (which has
  # always worked) handle the no-data case.
  observe({
    av <- availability()
    for (tab_id in names(av)) {
      avail <- av[[tab_id]]
      session$sendCustomMessage("set_disabled_options", list(
        input_id = paste0("group_", tab_id),
        # Wrap with I() so an empty vector serializes as `[]` (a bare empty
        # list becomes `{}`, which would break indexOf in the JS handler).
        disabled = I(as.character(names(avail)[!avail]))
      ))
    }
  })

  # Counts shown in the per-tab summary line. Trials and references come
  # from the survivor set (arms that actually contribute a row to the
  # endpoint); references are summed across pubs_by_study to count every
  # primary + secondary publication, not just unique study IDs.
  summarise_view <- function(survivors) {
    if (!nrow(survivors)) {
      return(list(n_trials = 0L, n_refs = 0L, n_patients = 0L))
    }
    refs <- unique(as.character(survivors$ref_id))
    n_refs <- sum(vapply(refs, function(rid) {
      p <- pubs_by_study[[rid]]
      if (is.null(p)) 0L else nrow(p)
    }, integer(1)))
    arms <- unique(survivors[, c("ref_id", "arm_no", "n_arm"), drop = FALSE])
    list(
      n_trials   = length(unique(survivors$trial)),
      n_refs     = n_refs,
      n_patients = sum(arms$n_arm, na.rm = TRUE)
    )
  }

  pluralise <- function(n, word) {
    sprintf("%s %s%s",
            formatC(n, format = "d", big.mark = ","),
            word,
            if (n == 1L) "" else "s")
  }

  # Build one renderDT per tab. Reactive picks the active endpoint group
  # from that tab's dropdown, queries the right table, formats it.
  for (tab_id in names(endpoint_groups)) local({
    this_tab    <- tab_id
    tab_cfg     <- endpoint_groups[[this_tab]]

    output[[paste0("summary_", this_tab)]] <- renderUI({
      gid <- input[[paste0("group_", this_tab)]]
      req(gid)
      grp <- tab_cfg$groups[[gid]]
      surv <- group_survivors[[this_tab]][[gid]](
        query_view(grp$table, filter_state()))
      s <- summarise_view(surv)
      div(class = "view-summary",
          tags$span(class = "view-summary-filter",
                    filter_text(filter_state())),
          tags$span(class = "view-summary-sep", " • "),
          sprintf("%s across %s • %s",
                  pluralise(s$n_trials,   "trial"),
                  pluralise(s$n_refs,     "publication"),
                  pluralise(s$n_patients, "patient")))
    })

    output[[paste0("tbl_", this_tab)]] <- renderDT({
      gid <- input[[paste0("group_", this_tab)]]
      req(gid)
      grp <- tab_cfg$groups[[gid]]
      df <- grp$fmt(query_view(grp$table, filter_state()))
      n_endpoint_cols <- ncol(df) - 3
      cap <- if (!is.null(grp$note))
        htmltools::tags$caption(
          style = "caption-side: top; text-align: left;
                   padding: 6px 4px; color: #5a6478; font-size: 12px;",
          grp$note
        )
      datatable(
        df,
        rownames  = FALSE,
        caption   = cap,
        filter    = "none",
        selection = "none",  # Row-highlight on click is sticky and noisy;
                              # we don't use selection for anything anyway.
        # Drop DT's default `stripe` class — banding is applied per-trial in
        # rowCallback below, not per-row.
        class     = "row-border hover order-column",
        escape    = -1,  # Only the Trial column (col 1) contains HTML (an
                        # <a href="https://doi.org/..."> link when a DOI is
                        # known); fmt_trial() escapes the visible trial name
                        # itself. Every other column stays escaped.
        options  = list(
          pageLength = 25,
          autoWidth  = FALSE,
          dom        = "tip",
          columnDefs = list(list(className = "dt-right",
                                 targets = 2:(2 + n_endpoint_cols))),
          # Banding by trial: walk the rendered page top-to-bottom and flip
          # a band index whenever the visible trial-cell text changes. CSS
          # below paints alternate trials with a subtle background.
          drawCallback = JS(
            "function() {",
            "  var rows = this.api().rows({page:'current'}).nodes();",
            "  var prev = null, band = 0;",
            "  rows.each(function(node) {",
            "    var t = node.cells[0].innerText;",
            "    if (t !== prev) { band = 1 - band; prev = t; }",
            "    node.classList.toggle('trial-band', band === 1);",
            "  });",
            "}"
          ),
          # On page change, scroll the surrounding .tab-content (our custom
          # scroll container) back to the top — DT only scrolls its own
          # internal viewport, which we don't use.
          initComplete = JS(
            "function() {",
            "  var api = this.api();",
            "  api.on('page.dt', function() {",
            "    var el = $(api.table().node()).closest('.tab-content')[0];",
            "    if (el) el.scrollTo({ top: 0 });",
            "  });",
            "}"
          )
        ),
        colnames = grp$colnames
      )
    })
  })

  # The active endpoint group's pre-built network (nodes + edges, with
  # master-graph x/y baked in). Subset switches don't relayout.
  current_network <- reactive({
    tid <- input$view %||% names(endpoint_groups)[1]
    gid <- input[[paste0("group_", tid)]]
    req(gid)
    group_networks[[tid]][[gid]]
  })

  output$nma <- renderVisNetwork({
    nw <- isolate(current_network())
    visNetwork(nw$nodes, nw$edges) |>
      visNodes(shape   = "dot",
               # Per-node `size` is set in build_network_data() (sqrt of
               # patient count, anchored to the global max). No `scaling`
               # block — that would re-rescale per render and break visual
               # comparability across endpoint subsets.
               font    = list(size = 44, face = "Helvetica",
                              strokeWidth = 6, strokeColor = "#ffffff"),
               color   = list(background = "#4C9AFF",
                              border     = "#1F4E8C",
                              highlight  = list(background = "#FF8A3D",
                                                border     = "#B5521A"),
                              hover      = list(background = "#7FB5FF",
                                                border     = "#1F4E8C")),
               borderWidth = 2,
               physics = FALSE) |>
      visEdges(smooth   = list(enabled = TRUE, type = "continuous"),
               # Floor the minimum width at 3px so single-trial edges stay
               # easily clickable; head-to-heads with many trials still
               # stand out at the top end.
               scaling  = list(min = 3, max = 14),
               color    = list(color = "rgba(80,80,80,0.35)",
                               highlight = "#FF8A3D",
                               hover     = "#1F4E8C")) |>
      visPhysics(enabled = FALSE) |>
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1,
                                         hover = TRUE,
                                         labelOnly = FALSE),
                 nodesIdSelection = FALSE) |>
      visInteraction(navigationButtons = FALSE, multiselect = FALSE,
                     tooltipDelay = 150, hover = TRUE,
                     zoomView = TRUE, dragView = TRUE) |>
      visEvents(
        selectNode   = "function(p){ Shiny.setInputValue('nma_node', p.nodes[0], {priority:'event'}); }",
        selectEdge   = "function(p){ if(p.nodes && p.nodes.length) return;
                                     Shiny.setInputValue('nma_edge', p.edges[0], {priority:'event'}); }",
        deselectNode = "function(p){ Shiny.setInputValue('nma_clear', Math.random(), {priority:'event'}); }",
        deselectEdge = "function(p){ Shiny.setInputValue('nma_clear', Math.random(), {priority:'event'}); }"
      )
  })

  # When the active endpoint group changes, swap nodes/edges in place via
  # the proxy so zoom and pan are preserved. If the currently-selected drug
  # or pair has no data in the new endpoint, clear the filter automatically
  # — the user confirmed this behaviour over leaving the selection pinned.
  observeEvent(current_network(), {
    nw <- current_network()
    remove_nodes <- setdiff(master_network$nodes$id, nw$nodes$id)
    remove_edges <- setdiff(master_network$edges$id, nw$edges$id)
    proxy <- visNetworkProxy("nma")
    if (length(remove_edges)) visRemoveEdges(proxy, id = remove_edges)
    if (length(remove_nodes)) visRemoveNodes(proxy, id = remove_nodes)
    visUpdateNodes(proxy, nodes = nw$nodes)
    visUpdateEdges(proxy, edges = nw$edges)
    # Recentre on the visible subset so a smaller network doesn't sit in a
    # corner of the canvas. Preserves zoom level on no-op updates.
    visFit(proxy, animation = list(duration = 250))

    s <- isolate(filter_state())
    if (!is.null(s)) {
      stale <- FALSE
      if (identical(s$kind, "node") && !(s$drug %in% nw$nodes$id)) {
        stale <- TRUE
      } else if (identical(s$kind, "edge")) {
        present <- nrow(nw$edges) > 0 && any(
          (nw$edges$from == s$from & nw$edges$to == s$to) |
          (nw$edges$from == s$to   & nw$edges$to == s$from)
        )
        if (!present) stale <- TRUE
      }
      if (stale) {
        filter_state(NULL)
        visUnselectAll(proxy)
      }
    }
  }, ignoreInit = TRUE)

  # Node click -> single-drug filter
  observeEvent(input$nma_node, {
    req(input$nma_node)
    filter_state(list(kind = "node", drug = input$nma_node))
  })

  # Edge click -> head-to-head pair filter
  observeEvent(input$nma_edge, {
    req(input$nma_edge)
    e <- edges_df[edges_df$id == input$nma_edge, , drop = FALSE]
    if (!nrow(e)) return()
    filter_state(list(kind = "edge", from = e$from[1], to = e$to[1]))
  })

  # Empty-canvas click -> clear
  observeEvent(input$nma_clear, { filter_state(NULL) })

  # Meta-analysis modal. Reads the precomputed ma_* tables and renders one
  # forest plot per outcome in the active endpoint group. Dispatch on
  # filter_state(): NULL -> drug vs placebo, edge -> pairwise, node -> single
  # arm proportion.
  ma_emphasis <- reactiveVal("RE")
  observeEvent(input$ma_emphasis, {
    req(input$ma_emphasis %in% c("FE", "RE"))
    ma_emphasis(input$ma_emphasis)
  })

  observeEvent(input$show_ma, {
    state  <- filter_state()
    tab_id <- input$view %||% names(endpoint_groups)[1]
    gid    <- input[[paste0("group_", tab_id)]]
    req(gid)

    if (!HAS_MA) {
      showModal(modalDialog(
        title = "Meta-analysis", easyClose = TRUE,
        footer = modalButton("Close"), class = "ma-modal",
        tags$p("No precomputed meta-analysis tables found. Run ",
               tags$code("Rscript app/meta_analyse.R"),
               " after ", tags$code("convert.R"),
               " to build them, then reload the app.")
      ))
      return(invisible(NULL))
    }

    outcomes <- ma_catalog[[tab_id]][[gid]]
    if (is.null(outcomes) || !length(outcomes)) {
      showModal(modalDialog(
        title = "Meta-analysis", easyClose = TRUE,
        footer = modalButton("Close"), class = "ma-modal",
        tags$p("Meta-analysis isn't configured for this endpoint group yet.")
      ))
      return(invisible(NULL))
    }

    # Resolve title + summary line.
    state_lbl <- if (is.null(state)) "All drugs vs placebo"
                 else if (identical(state$kind, "edge"))
                   sprintf("%s vs %s", state$from, state$to)
                 else if (identical(state$kind, "node"))
                   sprintf("Pooled response rate, %s", state$drug)
                 else "All drugs"
    group_lbl <- endpoint_groups[[tab_id]]$groups[[gid]]$label %||% gid

    # Build a plot block per outcome.
    plot_blocks <- lapply(outcomes, function(outc) {
      out_id <- paste0("ma_plot_", tab_id, "_", gid, "_", outc$code)
      inputs <- tryCatch(build_forest_inputs(state, tab_id, outc),
                        error = function(e) NULL)
      if (is.null(inputs)) {
        return(div(class = "ma-plot-block",
                   div(class = "ma-plot-title", outc$label),
                   div(class = "ma-empty",
                       "No meta-analysable data for this comparison.")))
      }
      # `effective_scale` lets the dispatcher override outc$scale -- e.g. a
      # node filter on a binary outcome plots single-arm proportions, not RR.
      eff_scale <- inputs$effective_scale %||% outc$scale
      stats <- if (!is.null(inputs$pooled)) {
        fe <- inputs$pooled[inputs$pooled$kind == "FE", ]
        re <- inputs$pooled[inputs$pooled$kind == "RE", ]
        sm_lbl <- switch(eff_scale, rr = "RR", md = "MD", prop = "Proportion")
        digits <- if (eff_scale == "prop") 3 else 2
        sprintf("FE: %s %s | RE: %s %s | %d studies",
                sm_lbl, fmt_ci(fe$est, fe$lo, fe$hi, digits),
                sm_lbl, fmt_ci(re$est, re$lo, re$hi, digits),
                nrow(inputs$rows))
      } else {
        sprintf("%d %s",
                nrow(inputs$rows),
                if (is.null(state)) "drugs" else "rows")
      }
      div(class = "ma-plot-block",
          div(class = "ma-plot-title", outc$label),
          div(class = "ma-plot-stats", stats),
          girafeOutput(out_id, height = "auto"))
    })

    showModal(modalDialog(
      title = tagList(icon("chart-column"), " Meta-analysis"),
      easyClose = TRUE, size = "l", class = "ma-modal",
      footer = tagList(
        span(class = "ma-footer-meta",
             if (!is.na(MA_BUILT_AT)) sprintf("Built: %s UTC", MA_BUILT_AT)
             else "Built: unknown"),
        span(class = "ma-toggle",
             "Emphasise:",
             radioButtons("ma_emphasis", label = NULL,
                          choices = c("Random effects" = "RE",
                                       "Fixed effect"   = "FE"),
                          selected = ma_emphasis(),
                          inline = TRUE)),
        modalButton("Close")
      ),
      div(
        div(class = "ma-summary",
            sprintf("%s | endpoint group: %s | random-effects pooling via REML",
                    state_lbl, group_lbl)),
        plot_blocks
      )
    ))

    # Register one renderer per plot. Reactive on ma_emphasis() so the diamond
    # colour swap is live.
    for (outc in outcomes) local({
      this_outc <- outc
      out_id    <- paste0("ma_plot_", tab_id, "_", gid, "_", this_outc$code)
      output[[out_id]] <- renderGirafe({
        inputs <- tryCatch(build_forest_inputs(state, tab_id, this_outc),
                          error = function(e) NULL)
        if (is.null(inputs)) return(NULL)
        forest_ggiraph(inputs$rows, inputs$pooled,
                       scale = inputs$effective_scale %||% this_outc$scale,
                       emphasis = ma_emphasis(),
                       title = NULL,
                       output_id = out_id)
      })
    })
  })

  observeEvent(input$show_about, {
    showModal(modalDialog(
      title = "About the Psoriasis RCT Explorer",
      easyClose = TRUE,
      size = "l",
      footer = modalButton("Close"),
      class = "about-modal",
      tags$div(
        tags$h4("Overview"),
        tags$p("This app can be used to explore the results from randomised 
               controlled trials in psoriasis. Data for psoriasis area and
               severity index (PASI), Dermatology Life Quality Index (DLQI),
               and safety are included. Use the network diagram on the left
               to choose the drug or the comparison you are interested in and
               the tables will update automatically."),
        tags$h4("Study selection"),
        tags$p("The studies reported in this app were identified as 
               randomised controlled trials for systemic treatments for
               moderate-to-severe psoriasis by a Cochrane living review ",
               tags$a("(Sbidian et al., 2025)", 
                      href = "https://doi.org/10.1002/14651858.CD011535.pub7",
                      target = "_blank", rel = "noopener", .noWS = "after"),
               ". "),
        tags$h4("Data extraction"),
        tags$p("Data were extracted from the identified studies by a single
               researcher. The data is stored in a sqlite database, which can
               be downloaded. For full details of the data extraction, see 
               the full protocol on ",
               tags$a("PROSPERO", 
                      href = "https://www.crd.york.ac.uk/PROSPERO/view/CRD420261306630",
                      target = "_blank", rel = "noopener", .noWS = "after"),
               "."
               ),
        tags$h4("Network diagram"),
        tags$p("The network diagram shows the drugs and comparisons available
               for the selected endpoints. Nodes are sizes according sample 
               size; edges are sizes according to number of trials."),
        tags$h4("Contact"),
        tags$p("If you have a question about the app, please raise an issue 
               on the ",
               tags$a("GitHub repository", 
                      href = "https://github.com/tomclemmet/psoriasis-rct-explorer",
                      target = "_blank", rel = "noopener", .noWS = "after"), 
               ".")
      )
    ))
  })

  output$download_db <- downloadHandler(
    filename    = function() "psoriasis-rcts.sqlite",
    contentType = "application/x-sqlite3",
    content     = function(file) file.copy(DB_PATH, file)
  )
}

shinyApp(ui, server)
