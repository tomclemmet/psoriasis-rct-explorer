suppressPackageStartupMessages({
  library(shiny)
  library(DBI)
  library(RSQLite)
  library(DT)
  library(visNetwork)
  library(jsonlite)
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
# "drug_vs_placebo"). Tooltips carry trial/drug name, n / N, effect and CI.
# Per-trial squares are sized by the trial's sample size; per-drug network
# rows are drawn uniform.
# ---------------------------------------------------------------------------
fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f (%.*f to %.*f)",
          digits, est, digits, lo, digits, hi)
}

# Plain-text tooltip for native SVG <title> elements. Multi-line via "\n";
# browsers honour newlines in <title> as a multi-line hover tooltip without
# any extra JS. `extra` is an optional named character vector of additional
# label / value pairs.
ma_tooltip <- function(label, est, lo, hi, extra = NULL, digits = 2) {
  parts <- c(
    label,
    sprintf("Effect: %s", fmt_ci(est, lo, hi, digits))
  )
  if (length(extra))
    parts <- c(parts, sprintf("%s: %s", names(extra), as.character(extra)))
  paste(parts, collapse = "\n")
}

# --- Inline-SVG forest plot ------------------------------------------------
# All rendering is a pure R string-builder: no ggplot, no ggiraph, no
# svglite. We compute pixel coordinates ourselves and emit <svg> with
# <line>, <rect>, <polygon>, <text> children. Per-row tooltips ride along
# as native SVG <title> child elements (every browser renders these as
# hover tooltips automatically).
#
# `rows` is a data.frame with columns
#   label, est, lo, hi, square_n, tooltip
# `square_n` is the trial's sample size used to size its square; NA on every
# row (e.g. network pooled rows) draws all squares uniform.
# Optional `pooled` rows are drawn as diamonds below; each has kind ("FE"
# or "RE"), est, lo, hi. `scale` is one of "rr", "md", "prop".

# Pick visible x-axis ticks for one of the three scales. Returns numeric
# tick values in data units; the renderer formats them and maps to pixels.
forest_ticks <- function(scale, xmin, xmax) {
  if (identical(scale, "rr")) {
    # Decade-anchored ticks; keep ones inside the visible window.
    candidates <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
    candidates[candidates >= xmin & candidates <= xmax]
  } else if (identical(scale, "prop")) {
    # Pick a tick interval that gives 4-6 ticks across the actual range.
    span <- xmax - xmin   # xmin is always 0
    step <- if      (span <= 0.10) 0.02
            else if (span <= 0.20) 0.05
            else if (span <= 0.50) 0.10
            else                   0.25
    seq(0, xmax, by = step)
  } else {
    # mean diff: ~5 pretty ticks within the data range.
    pretty(c(xmin, xmax), n = 5)
  }
}

# Format a tick label for the axis.
forest_tick_label <- function(scale, v) {
  if (identical(scale, "rr"))    sub("\\.0$", "", formatC(v, format = "g"))
  else if (identical(scale, "prop")) sprintf("%g", v)
  else                            formatC(v, format = "g")
}

# Map a data-unit x value to a pixel x position inside the plot area.
forest_xscale <- function(scale, xmin, xmax, plot_left, plot_w) {
  if (identical(scale, "rr")) {
    lmin <- log10(xmin); lmax <- log10(xmax)
    function(x) plot_left + (log10(pmax(x, xmin / 10)) - lmin) /
                            (lmax - lmin) * plot_w
  } else {
    function(x) plot_left + (x - xmin) / (xmax - xmin) * plot_w
  }
}

# Compute (xmin, xmax) so every CI fits with a little padding, and the
# reference line sits inside the visible range when it exists.
forest_xlimits <- function(scale, rows, pooled) {
  vals <- c(rows$lo, rows$hi)
  if (!is.null(pooled) && nrow(pooled)) vals <- c(vals, pooled$lo, pooled$hi)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(switch(scale, rr = c(0.1, 10), prop = c(0, 1),
                                   c(-1, 1)))
  if (identical(scale, "rr")) {
    # Pad in log space and snap outward to a decade-anchored tick so the
    # axis ends on a clean number. Always include 1 (the null) inside the
    # range so the reference line is visible.
    candidates <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
    lo <- min(vals, 1); hi <- max(vals, 1)
    # Add log-space padding before snapping so a CI that lands exactly on a
    # candidate tick doesn't have its tip touching the plot edge.
    lr  <- log10(hi) - log10(lo)
    pad <- max(lr * 0.05, 0.05)
    lo  <- 10 ^ (log10(lo) - pad)
    hi  <- 10 ^ (log10(hi) + pad)
    lo <- max(candidates[candidates <= lo], 0.01)
    hi <- min(candidates[candidates >= hi], 200)
    c(lo, hi)
  } else if (identical(scale, "prop")) {
    # Always start at 0. Upper bound: if all CIs fit below 0.5 use a tight
    # ceiling snapped to a clean fraction; otherwise fall back to 1.
    hi_data <- max(vals)
    if (hi_data <= 0.5) {
      # Snap upper bound outward to the nearest 0.05 increment, with a little
      # headroom so the highest CI tip doesn't sit at the edge.
      pad <- max(hi_data * 0.10, 0.02)
      hi_snap <- ceiling((hi_data + pad) / 0.05) * 0.05
      c(0, min(hi_snap, 0.5))
    } else {
      c(0, 1)
    }
  } else {
    pad <- (max(vals) - min(vals)) * 0.1
    # Linear (MD): also pull the null (0) inside the range when CIs lie on
    # one side of zero.
    c(min(min(vals), 0) - pad, max(max(vals), 0) + pad)
  }
}

# Diamond polygon: left point, top, right, bottom — half-height `hh` px.
forest_diamond_points <- function(xL, xC, xR, yC, hh) {
  sprintf("%g,%g %g,%g %g,%g %g,%g",
          xL, yC, xC, yC - hh, xR, yC, xC, yC + hh)
}

# Render one forest plot to a complete <svg> string. Width is fixed; height
# scales with the number of rows.
#
# `axis_label`  : single string drawn centred beneath the axis
#                 (e.g. "Risk ratio (95% CI, log scale)").
# `dir_left`    : optional small label drawn at the left end of the axis line
#                 (e.g. "favours placebo →" reversed).
# `dir_right`   : optional small label drawn at the right end.
forest_svg <- function(rows, pooled, scale = "rr", width = 880,
                       axis_label = NULL, dir_left = NULL, dir_right = NULL) {
  esc <- function(x) htmltools::htmlEscape(x, attribute = FALSE)
  esc_attr <- function(x) htmltools::htmlEscape(x, attribute = TRUE)
  if (!nrow(rows)) {
    return(sprintf('<div class="ma-empty">No meta-analysable data for this comparison.</div>'))
  }

  # Geometry.
  ROW_H        <- 22
  HEADER_PAD   <- 8
  POOLED_GAP   <- 14
  POOLED_H     <- 22
  AXIS_PAD     <- 28
  AXIS_LABEL_H <- if (length(axis_label) || length(dir_left) || length(dir_right)) 22 else 0
  BOTTOM_PAD   <- 6
  LEFT_MARGIN  <- 220       # row labels (a touch wider for "Pooled estimate" rows)
  RIGHT_MARGIN <- 175       # CI text only -- weight column dropped
  PLOT_LEFT    <- LEFT_MARGIN
  PLOT_W       <- width - LEFT_MARGIN - RIGHT_MARGIN

  n_trials  <- nrow(rows)
  n_pooled  <- if (!is.null(pooled)) nrow(pooled) else 0L

  body_h    <- n_trials * ROW_H +
    (if (n_pooled) POOLED_GAP + n_pooled * POOLED_H else 0)
  height    <- HEADER_PAD + body_h + AXIS_PAD + AXIS_LABEL_H + BOTTOM_PAD

  # Domain + tick positions.
  xlim  <- forest_xlimits(scale, rows, pooled)
  xmin  <- xlim[1]; xmax <- xlim[2]
  xfn   <- forest_xscale(scale, xmin, xmax, PLOT_LEFT, PLOT_W)
  ticks <- forest_ticks(scale, xmin, xmax)

  ref_line <- switch(scale, rr = 1, md = 0, prop = NA_real_)
  plot_bottom <- HEADER_PAD + body_h
  plot_top    <- HEADER_PAD - 4

  # CI digits: tighter for proportions because they tend to cluster.
  ci_digits <- if (identical(scale, "prop")) 3 else 2

  # Per-row square size from the trial's sample size; clamped to a readable
  # range. Rows without a sample size (e.g. network pooled rows) draw uniform.
  square_n <- if (!is.null(rows$square_n)) rows$square_n else rep(NA_real_, n_trials)
  max_n <- suppressWarnings(max(square_n, na.rm = TRUE))
  square_size <- if (!is.finite(max_n) || max_n <= 0) {
    rep(8, n_trials)
  } else {
    pmax(5, pmin(12, 5 + 7 * sqrt(square_n / max_n)))
  }
  square_size[is.na(square_size)] <- 8

  parts <- character(0)

  # Open the SVG element.
  parts[length(parts) + 1L] <- sprintf(
    '<svg class="ma-forest" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%" height="%dpx" preserveAspectRatio="xMinYMin meet" role="img">',
    width, as.integer(height), as.integer(height))

  # Plot-area background (subtle).
  parts[length(parts) + 1L] <- sprintf(
    '<rect class="ma-plot-bg" x="%g" y="%g" width="%g" height="%g"/>',
    PLOT_LEFT, plot_top, PLOT_W, plot_bottom - plot_top)

  # Reference line (RR=1, MD=0; nothing for proportions).
  if (!is.na(ref_line) && ref_line >= xmin && ref_line <= xmax) {
    rx <- xfn(ref_line)
    parts[length(parts) + 1L] <- sprintf(
      '<line class="ma-refline" x1="%g" y1="%g" x2="%g" y2="%g"/>',
      rx, plot_top, rx, plot_bottom)
  }

  # CI domain pixel boundaries (for clamping + off-scale arrows).
  X_MIN_PX <- PLOT_LEFT
  X_MAX_PX <- PLOT_LEFT + PLOT_W

  # Build a single trial / drug row. Squares clamp inside the axis;
  # CI lines that exceed bounds are drawn up to the edge and capped
  # with a triangular "off-scale" arrowhead.
  emit_data_row <- function(yc, est, lo, hi, sz, label_left, ci_text,
                            tooltip, klass = "ma-square-default") {
    xL_raw <- xfn(lo); xR_raw <- xfn(hi); xE_raw <- xfn(est)
    # NA-safe helpers: NA/NaN positions mean "no data", treated as in-range.
    finite_or <- function(x, fallback) if (is.finite(x)) x else fallback
    xL <- max(X_MIN_PX, finite_or(xL_raw, X_MIN_PX))
    xR <- min(X_MAX_PX, finite_or(xR_raw, X_MAX_PX))
    xE_raw_safe <- finite_or(xE_raw, (X_MIN_PX + X_MAX_PX) / 2)
    # Clamp the square's centre so its body stays inside the plot.
    xE_clamped <- min(max(xE_raw_safe, X_MIN_PX + sz / 2), X_MAX_PX - sz / 2)
    off_left  <- is.finite(xL_raw) && xL_raw < X_MIN_PX
    off_right <- is.finite(xR_raw) && xR_raw > X_MAX_PX
    has_data  <- is.finite(xE_raw)   # hide CI bar + square if estimate is NA
    paste0(
      '<g class="ma-row" data-tt="', esc_attr(tooltip), '">',
      sprintf('<text class="ma-rowlabel" x="%g" y="%g">%s</text>',
              LEFT_MARGIN - 10, yc + 4, esc(label_left)),
      # Hover catcher: full row, makes it easy to trigger the tooltip.
      sprintf('<rect class="ma-rowhit" x="%g" y="%g" width="%g" height="%g"/>',
              PLOT_LEFT, yc - ROW_H / 2 + 1, PLOT_W, ROW_H - 2),
      if (has_data && xL < xR)
        sprintf('<line class="ma-ci" x1="%g" y1="%g" x2="%g" y2="%g"/>',
                xL, yc, xR, yc) else "",
      if (off_left)
        sprintf('<polygon class="ma-ci-arrow" points="%g,%g %g,%g %g,%g"/>',
                X_MIN_PX, yc, X_MIN_PX + 6, yc - 4, X_MIN_PX + 6, yc + 4)
      else "",
      if (off_right)
        sprintf('<polygon class="ma-ci-arrow" points="%g,%g %g,%g %g,%g"/>',
                X_MAX_PX, yc, X_MAX_PX - 6, yc - 4, X_MAX_PX - 6, yc + 4)
      else "",
      if (has_data)
        sprintf('<rect class="ma-square %s" x="%g" y="%g" width="%g" height="%g"/>',
                klass, xE_clamped - sz / 2, yc - sz / 2, sz, sz)
      else "",
      sprintf('<text class="ma-citext" x="%g" y="%g">%s</text>',
              width - RIGHT_MARGIN + 10, yc + 4, esc(ci_text)),
      '</g>'
    )
  }

  # Optional per-row colour class (FE/RE in no-filter view); default = black.
  row_klass <- if (!is.null(rows$klass)) rows$klass
               else rep("ma-square-default", n_trials)

  for (i in seq_len(n_trials)) {
    yc <- HEADER_PAD + (i - 0.5) * ROW_H
    parts[length(parts) + 1L] <- emit_data_row(
      yc, rows$est[i], rows$lo[i], rows$hi[i],
      sz = square_size[i],
      label_left = rows$label[i],
      ci_text    = fmt_ci(rows$est[i], rows$lo[i], rows$hi[i], ci_digits),
      tooltip    = rows$tooltip[i] %||% rows$label[i],
      klass      = row_klass[i]
    )
  }

  # Pooled rows: diamond + label (kind + CI). Colour-coded by model only:
  # FE (pairwise or network) = mid-grey, RE = brand blue. Each row is labelled
  # on the left so pairwise vs network estimates are distinguished by text.
  if (n_pooled) {
    for (i in seq_len(n_pooled)) {
      yc <- HEADER_PAD + n_trials * ROW_H + POOLED_GAP + (i - 0.5) * POOLED_H

      xL_raw <- xfn(pooled$lo[i]); xR_raw <- xfn(pooled$hi[i])
      xE_raw <- xfn(pooled$est[i])
      xL <- max(X_MIN_PX, if (is.finite(xL_raw)) xL_raw else X_MIN_PX)
      xR <- min(X_MAX_PX, if (is.finite(xR_raw)) xR_raw else X_MAX_PX)
      xE <- if (is.finite(xE_raw)) min(max(xE_raw, X_MIN_PX), X_MAX_PX)
            else (X_MIN_PX + X_MAX_PX) / 2
      kind  <- pooled$kind[i]
      klass <- switch(kind,
        "FE"     = "ma-pooled-fe",
        "RE"     = "ma-pooled-re",
        "NMA-FE" = "ma-pooled-fe",
        "NMA-RE" = "ma-pooled-re",
        "ma-pooled-fe")
      kind_lbl <- switch(kind,
        "FE"     = "Pooled estimate (FE)",
        "RE"     = "Pooled estimate (RE)",
        "NMA-FE" = "Network estimate (FE)",
        "NMA-RE" = "Network estimate (RE)",
        kind)
      pts    <- forest_diamond_points(xL, xE, xR, yc, 7)
      ci_str <- fmt_ci(pooled$est[i], pooled$lo[i], pooled$hi[i], ci_digits)
      tt_hdr <- switch(kind,
        "FE"     = "Common-effect (FE) pooled estimate",
        "RE"     = "Random-effects (RE, REML) pooled estimate",
        "NMA-FE" = "Network MA common-effect (FE) estimate",
        "NMA-RE" = "Network MA random-effects (RE, REML) estimate",
        kind)
      tt_extra <- if (!is.null(pooled$n_direct) && !is.na(pooled$n_direct[i]))
        sprintf("%s\nDirect studies: %d", ci_str, pooled$n_direct[i])
      else
        ci_str
      tt <- sprintf("%s\n%s", tt_hdr, tt_extra)
      parts[length(parts) + 1L] <- paste0(
        sprintf('<g class="ma-pooled %s" data-tt="%s">', klass, esc_attr(tt)),
        sprintf('<text class="ma-rowlabel ma-pooled-label" x="%g" y="%g">%s</text>',
                LEFT_MARGIN - 10, yc + 4, kind_lbl),
        sprintf('<rect class="ma-rowhit" x="%g" y="%g" width="%g" height="%g"/>',
                PLOT_LEFT, yc - POOLED_H / 2 + 1, PLOT_W, POOLED_H - 2),
        sprintf('<polygon class="ma-diamond" points="%s"/>', pts),
        sprintf('<text class="ma-citext" x="%g" y="%g">%s</text>',
                width - RIGHT_MARGIN + 10, yc + 4, esc(ci_str)),
        '</g>'
      )
    }
  }

  # X-axis: a horizontal line + tick marks + labels.
  axis_y <- plot_bottom + 6
  parts[length(parts) + 1L] <- sprintf(
    '<line class="ma-axis" x1="%g" y1="%g" x2="%g" y2="%g"/>',
    PLOT_LEFT, axis_y, PLOT_LEFT + PLOT_W, axis_y)
  for (t in ticks) {
    tx <- xfn(t)
    parts[length(parts) + 1L] <- sprintf(
      '<line class="ma-tick" x1="%g" y1="%g" x2="%g" y2="%g"/>',
      tx, axis_y, tx, axis_y + 4)
    parts[length(parts) + 1L] <- sprintf(
      '<text class="ma-ticklabel" x="%g" y="%g">%s</text>',
      tx, axis_y + 16, esc(forest_tick_label(scale, t)))
  }

  # Axis label / directional cues, drawn below the tick labels.
  if (AXIS_LABEL_H > 0) {
    label_y <- axis_y + 30
    if (length(dir_left) && nzchar(dir_left)) {
      parts[length(parts) + 1L] <- sprintf(
        '<text class="ma-dir ma-dir-left" x="%g" y="%g">%s</text>',
        PLOT_LEFT, label_y, esc(dir_left))
    }
    if (length(dir_right) && nzchar(dir_right)) {
      parts[length(parts) + 1L] <- sprintf(
        '<text class="ma-dir ma-dir-right" x="%g" y="%g">%s</text>',
        PLOT_LEFT + PLOT_W, label_y, esc(dir_right))
    }
    if (length(axis_label) && nzchar(axis_label)) {
      parts[length(parts) + 1L] <- sprintf(
        '<text class="ma-axislabel" x="%g" y="%g">%s</text>',
        PLOT_LEFT + PLOT_W / 2, label_y, esc(axis_label))
    }
  }

  parts[length(parts) + 1L] <- "</svg>"
  paste(parts, collapse = "")
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

# Network meta-analysis: one row per drug-vs-placebo NMA estimate plus
# network-level heterogeneity / inconsistency stats. Returns a list with
#   $status   "ok" | "sparse" | "missing"  (missing = network row not found)
#   $estimates  data frame of drug-vs-placebo rows (empty when sparse/missing)
#   $summary    one-row data frame of network-level stats
fetch_nma_vs_placebo <- function(endpoint_group, outcome_code) {
  n <- read_db(
    "SELECT * FROM ma_nma WHERE endpoint_group = ? AND outcome_code = ?",
    params = list(endpoint_group, outcome_code))
  if (!nrow(n)) return(list(status = "missing"))
  if (!identical(n$status[1], "ok"))
    return(list(status = n$status[1], summary = n[1, , drop = FALSE]))
  est <- read_db(
    "SELECT e.* FROM ma_nma_estimates e
      WHERE e.network_id = ? AND e.drug_b = 'Placebo' AND e.drug_a <> 'Placebo'
      ORDER BY e.drug_a",
    params = list(n$network_id[1]))
  list(status = "ok", summary = n[1, , drop = FALSE], estimates = est)
}

# NMA estimate for a specific drug pair from one network. Returns NULL when
# the pair isn't in the network or the network is absent/sparse. `flipped`
# is TRUE when the stored row is drug_b vs drug_a; callers must negate log
# estimates and swap lo/hi before use.
fetch_nma_for_pair <- function(drug_a, drug_b, endpoint_group, outcome_code) {
  n <- read_db(
    "SELECT * FROM ma_nma WHERE endpoint_group = ? AND outcome_code = ?",
    params = list(endpoint_group, outcome_code))
  if (!nrow(n) || !identical(n$status[1], "ok")) return(NULL)
  nid <- n$network_id[1]
  for (orient in list(c(drug_a, drug_b), c(drug_b, drug_a))) {
    e <- read_db(
      "SELECT * FROM ma_nma_estimates WHERE network_id = ? AND drug_a = ? AND drug_b = ?",
      params = list(nid, orient[1], orient[2]))
    if (nrow(e)) {
      return(list(est     = e[1, , drop = FALSE],
                  flipped = !identical(orient[1], drug_a),
                  summary = n[1, , drop = FALSE]))
    }
  }
  NULL
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
build_forest_inputs <- function(state, tab_id, outcome, ma_kind = "RE") {
  endpoint_group <- ma_endpoint_group(tab_id, outcome)
  digits <- if (identical(outcome$scale, "md")) 2 else 2

  # Direction semantics. For benefit endpoints (PASI/DLQI responders),
  # higher = better, so RR > 1 favours drug_a. For harm endpoints (safety),
  # lower = better, so RR < 1 favours drug_a — directions flip.
  is_harm <- identical(tab_id, "safety")

  # Axis text presets per scale, refined below for specific filter contexts.
  axis_label <- switch(outcome$scale,
                       rr   = "Risk ratio (95% CI, log scale)",
                       md   = "Mean difference (95% CI)",
                       prop = "Proportion (95% CI)",
                       "Effect (95% CI)")

  if (is.null(state)) {
    # No filter -> network meta-analysis vs Placebo. One row per drug, single
    # estimate (FE or RE, controlled by `ma_kind`).
    nma <- fetch_nma_vs_placebo(endpoint_group, outcome$code)
    if (identical(nma$status, "missing")) return(NULL)
    if (identical(nma$status, "sparse")) {
      return(list(empty_reason = paste0(
        "Network too sparse for meta-analysis on this outcome ",
        "(no connected network with Placebo)."),
        nma_summary = nma$summary))
    }
    df <- nma$estimates
    if (!nrow(df)) return(NULL)
    bt <- function(x) if (identical(outcome$scale, "rr")) exp(x) else x
    use_re <- identical(ma_kind, "RE")
    est <- bt(if (use_re) df$te_re else df$te_fe)
    lo  <- bt(if (use_re) df$lo_re else df$lo_fe)
    hi  <- bt(if (use_re) df$hi_re else df$hi_fe)
    rows <- data.frame(
      label     = df$drug_a,
      est       = est, lo = lo, hi = hi,
      square_n  = NA_real_,
      klass     = if (use_re) "ma-square-re" else "ma-square-fe",
      stringsAsFactors = FALSE
    )
    kind_lbl <- if (use_re) "Random effects" else "Common effect"
    ns       <- nma$summary
    i2_str   <- if (!is.na(ns$i2)) sprintf("%.1f%%", 100 * ns$i2) else "n/a"
    rows$tooltip <- vapply(seq_len(nrow(rows)), function(i) {
      ma_tooltip(
        sprintf("%s vs Placebo — NMA pooled (%s)",
                df$drug_a[i], if (use_re) "RE" else "FE"),
        rows$est[i], rows$lo[i], rows$hi[i],
        extra = c(
          "Direct studies"  = as.character(df$n_direct[i]),
          "Network studies" = as.character(ns$n_studies),
          "Network I-squared" = i2_str))
    }, character(1))
    if (identical(outcome$scale, "rr")) {
      if (is_harm) {
        dir_left  <- "← favours drug"
        dir_right <- "favours placebo →"
      } else {
        dir_left  <- "← favours placebo"
        dir_right <- "favours drug →"
      }
    } else if (identical(outcome$scale, "md")) {
      # MD: PASI/DLQI lower is better, so MD < 0 favours drug.
      dir_left  <- "← favours drug"
      dir_right <- "favours placebo →"
    } else {
      dir_left <- dir_right <- NULL
    }
    return(list(rows = rows, pooled = NULL,
                axis_label = axis_label,
                dir_left = dir_left, dir_right = dir_right,
                nma_summary = ns, ma_kind = ma_kind))
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
      square_n = t$n_a + t$n_b,
      stringsAsFactors = FALSE
    )
    pooled <- data.frame(
      kind    = c("FE", "RE"),
      est     = c(s$te_fe, s$te_re),
      lo      = c(s$lo_fe, s$lo_re),
      hi      = c(s$hi_fe, s$hi_re),
      n_direct = NA_integer_,
      stringsAsFactors = FALSE
    )
    # Append NMA estimate for this pair (both FE and RE) when available.
    # Orient the lookup against the pairwise display direction (drug_a vs
    # drug_b) so the NMA estimate matches the per-trial rows, not state$from/to.
    nma_pair <- tryCatch(
      fetch_nma_for_pair(s$drug_a, s$drug_b, endpoint_group, outcome$code),
      error = function(e) NULL)
    if (!is.null(nma_pair)) {
      e <- nma_pair$est
      if (nma_pair$flipped) {
        # Negate log estimate and swap lo/hi to restore the drug_a→drug_b
        # direction used by the pairwise rows.
        nma_pooled <- data.frame(
          kind    = c("NMA-FE", "NMA-RE"),
          est     = c(-e$te_fe, -e$te_re),
          lo      = c(-e$hi_fe, -e$hi_re),
          hi      = c(-e$lo_fe, -e$lo_re),
          n_direct = as.integer(e$n_direct),
          stringsAsFactors = FALSE
        )
      } else {
        nma_pooled <- data.frame(
          kind    = c("NMA-FE", "NMA-RE"),
          est     = c(e$te_fe, e$te_re),
          lo      = c(e$lo_fe, e$lo_re),
          hi      = c(e$hi_fe, e$hi_re),
          n_direct = as.integer(e$n_direct),
          stringsAsFactors = FALSE
        )
      }
      pooled <- rbind(pooled, nma_pooled)
    }
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
                 extra = extra)
    }, character(1))
    dir_left  <- dir_right <- NULL
    if (identical(outcome$scale, "rr")) {
      if (is_harm) {
        dir_left  <- sprintf("← favours %s", s$drug_a)
        dir_right <- sprintf("favours %s →", s$drug_b)
      } else {
        dir_left  <- sprintf("← favours %s", s$drug_b)
        dir_right <- sprintf("favours %s →", s$drug_a)
      }
    } else if (identical(outcome$scale, "md")) {
      # PASI / DLQI absolute scores: lower is better, so a negative MD
      # (drug_a − drug_b < 0) favours drug_a.
      dir_left  <- sprintf("← favours %s", s$drug_a)
      dir_right <- sprintf("favours %s →", s$drug_b)
    }
    return(list(rows = rows, pooled = pooled,
                comparison = sprintf("%s vs %s", s$drug_a, s$drug_b),
                axis_label = axis_label,
                dir_left = dir_left, dir_right = dir_right))
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
      square_n = t$n,
      stringsAsFactors = FALSE
    )
    event_lbl <- if (is_harm) "Events" else "Responders"
    rows$tooltip <- vapply(seq_len(nrow(rows)), function(i)
      ma_tooltip(rows$label[i], rows$est[i], rows$lo[i], rows$hi[i],
                 extra = c(setNames(sprintf("%d / %d", t$k[i], t$n[i]),
                                    event_lbl)),
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
    prop_lbl <- if (is_harm) "Event rate" else "Response rate"
    return(list(rows = rows, pooled = pooled, drug = state$drug,
                effective_scale = "prop",
                axis_label = sprintf("%s, %s (95%% CI)", prop_lbl, state$drug),
                dir_left = NULL, dir_right = NULL))
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

    // Fast-popup tooltip for forest plot rows. Triggered on mouseover of any
    // element carrying data-tt; positioned next to the cursor with a tiny
    // offset and dismissed on mouseleave. No delay (the browser's native
    // <title> attribute has a ~500 ms hover delay we want to avoid).
    (function() {
      var tt = null;
      function ensure() {
        if (tt) return tt;
        tt = document.createElement('div');
        tt.className = 'ma-tt-popover';
        tt.style.display = 'none';
        document.body.appendChild(tt);
        return tt;
      }
      function hide() { if (tt) tt.style.display = 'none'; }
      function show(el, ev) {
        var t = ensure();
        t.textContent = el.getAttribute('data-tt') || '';
        t.style.display = 'block';
        var pw = t.offsetWidth, ph = t.offsetHeight;
        var x = ev.clientX + 14, y = ev.clientY + 14;
        if (x + pw > window.innerWidth - 8)  x = ev.clientX - pw - 14;
        if (y + ph > window.innerHeight - 8) y = ev.clientY - ph - 14;
        t.style.left = (window.scrollX + Math.max(4, x)) + 'px';
        t.style.top  = (window.scrollY + Math.max(4, y)) + 'px';
      }
      document.addEventListener('mouseover', function(ev) {
        var el = ev.target.closest && ev.target.closest('[data-tt]');
        if (el) show(el, ev);
      });
      document.addEventListener('mousemove', function(ev) {
        var el = ev.target.closest && ev.target.closest('[data-tt]');
        if (el && tt && tt.style.display !== 'none') show(el, ev);
        else if (!el) hide();
      });
      document.addEventListener('mouseout', function(ev) {
        if (!ev.relatedTarget ||
            !(ev.relatedTarget.closest && ev.relatedTarget.closest('[data-tt]')))
          hide();
      });
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
      background: #ffffff; padding: 9px 4px;
      min-height: 52px; box-sizing: border-box;
      font-size: 16px; font-weight: 500; color: #1a1f2c;
      border-bottom: 1px solid #e3e6ea;
      display: flex; align-items: center; gap: 12px;
    }
    .view-summary-text { flex: 1 1 auto; min-width: 0; }
    .view-summary .ma-btn {
      flex: 0 0 auto; align-self: center;
      padding: 6px 14px; font-size: 13px; font-weight: 600;
      background: #1F4E8C; border-color: #1F4E8C; color: #ffffff;
      box-shadow: 0 1px 2px rgba(0,0,0,0.08);
    }
    .view-summary .ma-btn:hover,
    .view-summary .ma-btn:focus {
      background: #163a6b; border-color: #163a6b; color: #ffffff;
    }
    .view-summary .ma-btn .fa { margin-right: 4px; }
    .row.split .col-table .tab-content table.dataTable thead th {
      position: sticky; top: 98px; background: #ffffff; z-index: 2;
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
    /* FE/RE toggle row (no-filter NMA view only). */
    .ma-modal .ma-toggle-row {
      display: flex; align-items: center; gap: 12px;
      margin: 4px 0 16px 0; padding: 7px 12px;
      background: #f6f7f9; border: 1px solid #e3e6ea; border-radius: 6px;
    }
    .ma-modal .ma-toggle-label {
      font-size: 12px; font-weight: 600; color: #5a6478;
      white-space: nowrap; flex: 0 0 auto; line-height: 1;
    }
    /* Strip Bootstrap's default form-group bottom margin so flex alignment works */
    .ma-modal .ma-toggle-row .shiny-input-container {
      margin-bottom: 0; flex: 0 0 auto;
    }
    .ma-modal .ma-toggle-row .shiny-options-group {
      display: flex; flex-direction: row; align-items: center;
      gap: 0 16px; margin: 0;
    }
    .ma-modal .ma-toggle-row .radio-inline {
      font-size: 13px; color: #1a1f2c;
      margin: 0; padding-left: 0; line-height: 1;
    }
    .ma-modal .ma-toggle-row .radio-inline input[type=radio] {
      margin: 0 5px 0 0; vertical-align: middle;
      position: relative; top: -1px;
    }
    .ma-modal .ma-footer-meta { font-size: 12px; color: #888;
                                margin-right: auto; }
    /* Inline-SVG forest plot styling. All forest plots share these rules;
       the SVG content itself is produced by forest_svg() (no ggplot). */
    .ma-forest { display: block; max-width: 100%; height: auto;
                 font-family: inherit; }
    .ma-forest .ma-plot-bg   { fill: #fafbfc; stroke: none; }
    .ma-forest .ma-refline   { stroke: #aab1bd; stroke-width: 1;
                                stroke-dasharray: 4 4; }
    .ma-forest .ma-axis      { stroke: #5a6478; stroke-width: 1; }
    .ma-forest .ma-tick      { stroke: #5a6478; stroke-width: 1; }
    .ma-forest .ma-ticklabel { fill: #5a6478; font-size: 11px;
                                text-anchor: middle; }
    .ma-forest .ma-rowlabel  { fill: #1a1f2c; font-size: 12px;
                                text-anchor: end; }
    .ma-forest .ma-ci        { stroke: #5a6478; stroke-width: 1.2; }
    .ma-forest .ma-ci-arrow  { fill: #5a6478; stroke: none; }
    .ma-forest .ma-square    { stroke: none; }
    .ma-forest .ma-square-default { fill: #1a1f2c; }
    .ma-forest .ma-square-fe { fill: #7f8fa6; }
    .ma-forest .ma-square-re { fill: #1F4E8C; }
    .ma-forest .ma-citext    { fill: #1a1f2c; font-size: 11px;
                                text-anchor: start; }
    .ma-forest .ma-pooled-fe .ma-diamond { fill: #7f8fa6; stroke: #5a6478;
                                            stroke-width: 1; }
    .ma-forest .ma-pooled-re .ma-diamond { fill: #1F4E8C; stroke: #14376a;
                                            stroke-width: 1; }
    .ma-forest .ma-pooled-label { font-weight: 600; }
    .ma-forest .ma-axislabel { fill: #1a1f2c; font-size: 12px;
                                text-anchor: middle; font-weight: 500; }
    .ma-forest .ma-dir       { fill: #5a6478; font-size: 11px;
                                font-style: italic; }
    .ma-forest .ma-dir-left  { text-anchor: start; }
    .ma-forest .ma-dir-right { text-anchor: end; }
    /* Invisible row hit-area: gives the cursor a wide target so the tooltip
       fires almost anywhere on the row, not just over the 8 px square. */
    .ma-forest .ma-rowhit { fill: transparent; pointer-events: all;
                            cursor: default; }
    .ma-forest .ma-row:hover .ma-square,
    .ma-forest .ma-row:hover .ma-ci   { stroke: #c0392b; }
    .ma-forest .ma-row:hover .ma-square-default { fill: #c0392b; }
    .ma-forest .ma-pooled:hover .ma-diamond { stroke: #c0392b; stroke-width: 1.5; }
    /* Fast-popup tooltip (JS-driven; replaces native <title>). */
    .ma-tt-popover {
      position: absolute; z-index: 1100; pointer-events: none;
      background: #1a1f2c; color: #f6f7f9; font-size: 12px;
      line-height: 1.4; padding: 6px 9px; border-radius: 4px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.18); white-space: pre-line;
      max-width: 320px;
    }
  "))),
  div(class = "title-bar",
      titlePanel("Psoriasis RCT Explorer"),
      div(class = "title-actions",
          downloadButton("download_db", "Download SQLite",
                         class = "btn btn-default"),
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
          tags$span(class = "view-summary-text",
              tags$span(class = "view-summary-filter",
                        filter_text(filter_state())),
              tags$span(class = "view-summary-sep", " • "),
              sprintf("%s across %s • %s",
                      pluralise(s$n_trials,   "trial"),
                      pluralise(s$n_refs,     "publication"),
                      pluralise(s$n_patients, "patient"))),
          actionButton(paste0("show_ma_", this_tab), "Meta-analysis",
                       icon = icon("chart-column"),
                       class = "btn btn-sm ma-btn"))
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

  # Meta-analysis modal. All forest plots are rendered as plain inline SVG
  # at modal-open time -- no Shiny output bindings, no ggplot, no ggiraph.
  # No-filter view shows a network meta-analysis vs Placebo with a
  # user-controlled FE/RE switch; edge/node views still draw both FE and RE
  # diamonds together (no switch).
  ma_ctx <- reactiveValues(active = FALSE, tab_id = NULL, state = NULL,
                           gid = NULL, outcomes = NULL, state_lbl = NULL,
                           group_lbl = NULL)
  ma_kind_state <- reactiveValues()  # per-tab last-selected ("RE"/"FE")

  # Build a single plot block for one outcome under the active state.
  build_plot_block <- function(state, tab_id, outc, ma_kind) {
    inputs <- tryCatch(build_forest_inputs(state, tab_id, outc, ma_kind),
                       error = function(e) NULL)
    if (is.null(inputs)) {
      return(div(class = "ma-plot-block",
                 div(class = "ma-plot-title", outc$label),
                 div(class = "ma-empty",
                     "No meta-analysable data for this comparison.")))
    }
    # Sparse-network sentinel from the no-filter NMA branch.
    if (!is.null(inputs$empty_reason)) {
      return(div(class = "ma-plot-block",
                 div(class = "ma-plot-title", outc$label),
                 div(class = "ma-empty", inputs$empty_reason)))
    }
    eff_scale <- inputs$effective_scale %||% outc$scale
    stats <- if (!is.null(inputs$pooled)) {
      fe <- inputs$pooled[inputs$pooled$kind == "FE", ]
      re <- inputs$pooled[inputs$pooled$kind == "RE", ]
      sm_lbl <- switch(eff_scale, rr = "RR", md = "MD", prop = "Proportion")
      digits <- if (eff_scale == "prop") 3 else 2
      sprintf("Pooled %s — FE: %s • RE: %s • %d studies",
              sm_lbl,
              fmt_ci(fe$est, fe$lo, fe$hi, digits),
              fmt_ci(re$est, re$lo, re$hi, digits),
              nrow(inputs$rows))
    } else if (!is.null(inputs$nma_summary)) {
      ns <- inputs$nma_summary
      kind_lbl <- if (identical(ma_kind, "FE")) "fixed effects"
                  else "random effects (REML)"
      i2_str <- if (!is.na(ns$i2)) sprintf("%.1f%%", 100 * ns$i2) else "n/a"
      pinc <- if (!is.na(ns$p_inc)) sprintf("%.2f", ns$p_inc) else "n/a"
      sprintf("Network MA, %s — %d studies, %d treatments • I² = %s • inconsistency p = %s",
              kind_lbl, ns$n_studies, ns$n_treatments, i2_str, pinc)
    } else {
      sprintf("%d %s",
              nrow(inputs$rows),
              if (is.null(state)) "drugs" else "rows")
    }
    svg_html <- forest_svg(inputs$rows, inputs$pooled, scale = eff_scale,
                           axis_label = inputs$axis_label,
                           dir_left   = inputs$dir_left,
                           dir_right  = inputs$dir_right)
    div(class = "ma-plot-block",
        div(class = "ma-plot-title", outc$label),
        div(class = "ma-plot-stats", stats),
        HTML(svg_html))
  }

  # Reactively rendered plot area. Driven by ma_ctx (set on modal open) and
  # input$ma_kind (toggle). Edge/node views ignore ma_kind.
  output$ma_plot_area <- renderUI({
    req(ma_ctx$active, ma_ctx$tab_id, ma_ctx$gid, ma_ctx$outcomes)
    tab_id   <- ma_ctx$tab_id
    state    <- ma_ctx$state
    outcomes <- ma_ctx$outcomes
    ma_kind  <- input$ma_kind %||% (ma_kind_state[[tab_id]] %||% "RE")
    lapply(outcomes, function(outc) build_plot_block(state, tab_id, outc, ma_kind))
  })

  # Persist toggle selection per tab across modal opens.
  observeEvent(input$ma_kind, {
    if (!is.null(ma_ctx$tab_id) && nzchar(input$ma_kind))
      ma_kind_state[[ma_ctx$tab_id]] <- input$ma_kind
  })

  open_ma_modal <- function(tab_id) {
    state  <- filter_state()
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

    state_lbl <- if (is.null(state)) "Network meta-analysis vs Placebo"
                 else if (identical(state$kind, "edge"))
                   sprintf("%s vs %s", state$from, state$to)
                 else if (identical(state$kind, "node"))
                   sprintf("Pooled response rate, %s", state$drug)
                 else "All drugs"
    group_lbl <- endpoint_groups[[tab_id]]$groups[[gid]]$label %||% gid

    # Stash context for the reactive plot area to consume.
    ma_ctx$tab_id    <- tab_id
    ma_ctx$state     <- state
    ma_ctx$gid       <- gid
    ma_ctx$outcomes  <- outcomes
    ma_ctx$state_lbl <- state_lbl
    ma_ctx$group_lbl <- group_lbl
    ma_ctx$active    <- TRUE

    summary_text <- if (is.null(state))
      sprintf("%s | endpoint group: %s | toggle FE/RE below",
              state_lbl, group_lbl)
    else
      sprintf("%s | endpoint group: %s | both common (FE) and random (RE, REML) pooled estimates shown",
              state_lbl, group_lbl)

    toggle_ui <- if (is.null(state)) {
      cur_kind <- ma_kind_state[[tab_id]] %||% "RE"
      div(class = "ma-toggle-row",
          tags$span(class = "ma-toggle-label", "Network model:"),
          radioButtons("ma_kind", label = NULL,
                       choices = c("Random effects" = "RE",
                                   "Fixed effects"  = "FE"),
                       selected = cur_kind, inline = TRUE))
    } else NULL

    showModal(modalDialog(
      title = tagList(icon("chart-column"), " Meta-analysis"),
      easyClose = TRUE, size = "l", class = "ma-modal",
      footer = tagList(
        span(class = "ma-footer-meta",
             if (!is.na(MA_BUILT_AT)) sprintf("Built: %s UTC", MA_BUILT_AT)
             else "Built: unknown"),
        modalButton("Close")
      ),
      div(
        div(class = "ma-summary", summary_text),
        toggle_ui,
        uiOutput("ma_plot_area")
      )
    ))
  }


  # One observer per tab — the button lives inside each tab's summary
  # output (so it's visually associated with the active filter context).
  for (tab_id in names(endpoint_groups)) local({
    this_tab <- tab_id
    observeEvent(input[[paste0("show_ma_", this_tab)]], {
      open_ma_modal(this_tab)
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
