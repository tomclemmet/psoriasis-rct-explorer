suppressPackageStartupMessages({
  library(shiny)
  library(DBI)
  library(RSQLite)
  library(DT)
  library(visNetwork)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

DERIVE_CHANGE_CORR <- 0.5  # assumed r(baseline, follow-up) for SD derivation

DB_PATH <- file.path(dirname(sys.frame(1)$ofile %||% "."), "psoriasis-rcts.sqlite")
if (!file.exists(DB_PATH)) DB_PATH <- "app/psoriasis-rcts.sqlite"
if (!file.exists(DB_PATH)) stop("psoriasis-rcts.sqlite not found - run convert.R first.")

source("db.R",     local = TRUE)
source("forest.R", local = TRUE)

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

.study_names <- read_db("SELECT study_id, trial FROM studies")
trial_name <- setNames(.study_names$trial, as.character(.study_names$study_id))

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

# Pre-rendered citation HTML per study_id, used by the trial detail modal.
cites_html_by_study <- vapply(
  pubs_by_study, fmt_citations, character(1)
)
cites_html_by_study <- cites_html_by_study[nzchar(cites_html_by_study)]
studies_with_cites <- names(cites_html_by_study)

# Render trial text as a clickable link. Clicking sends the ref-id to Shiny
# which opens a detail modal. Trials with no publications get a plain span.
# Caller must pass `escape = FALSE` for the Trial column.
fmt_trial <- function(trial, ref_id) {
  labels <- htmltools::htmlEscape(trial)
  ids    <- as.character(ref_id)
  has    <- ids %in% studies_with_cites
  out    <- sprintf("<span>%s</span>", labels)
  out[has] <- sprintf(
    '<a href="javascript:void(0)" class="trial-pop" data-ref-id="%s">%s</a>',
    ids[has], labels[has]
  )
  out
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

# Build a two-row arm header for the trial modal's tables: drug name on top
# (spanning consecutive arms sharing that drug), dose + sample size below.
# `dose_label` is the per-arm second-row text, e.g. "45 mg (n = 209)".
build_arm_header <- function(drug, dose_label) {
  esc  <- htmltools::htmlEscape
  grp  <- rle(drug)
  row1 <- paste0(
    sprintf('<th colspan="%d">%s</th>', grp$lengths, esc(grp$values)),
    collapse = ""
  )
  row2 <- paste0(sprintf("<th>%s</th>", esc(dose_label)), collapse = "")
  paste0(
    "<tr><th rowspan=\"2\"></th>", row1, "</tr>",
    "<tr>", row2, "</tr>"
  )
}

# Fill in the change column when both baseline and follow-up are reported
# but the study didn't report change directly. Returns the (possibly
# updated) change vector and a logical flag per row marking derived cells.
# Only the mean is derived: the SD of the difference depends on the
# within-arm baseline/follow-up correlation, which is rarely reported, so
# SDs of derived cells are left blank rather than guessed.
derive_change <- function(baseline, follow, change,
                          baseline_sd = NULL, follow_sd = NULL) {
  cd  <- rep(FALSE, length(change))
  csd <- rep(NA_real_, length(change))
  m   <- is.na(change) & !is.na(baseline) & !is.na(follow)
  change[m] <- follow[m] - baseline[m]
  cd[m]     <- TRUE
  if (!is.null(baseline_sd) && !is.null(follow_sd)) {
    can      <- m & !is.na(baseline_sd) & !is.na(follow_sd)
    csd[can] <- sqrt(
      follow_sd[can]^2 + baseline_sd[can]^2 -
        2 * DERIVE_CHANGE_CORR * follow_sd[can] * baseline_sd[can]
    )
  }
  list(change = change, change_derived = cd, change_sd_derived = csd)
}

# fmt_mean_sd, but derived cells are rendered without an SD (we don't have
# one) and tagged with a trailing asterisk so the table caption can explain
# them. Non-derived cells render exactly as fmt_mean_sd would.
fmt_mean_sd_marked <- function(mean, sd, derived, derived_sd = NULL, digits = 1) {
  sd_eff <- sd
  if (!is.null(derived_sd)) {
    fill <- derived & is.na(sd_eff) & !is.na(derived_sd)
    sd_eff[fill] <- derived_sd[fill]
  }
  out  <- fmt_mean_sd(mean, sd_eff, digits = digits)
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
  # Baseline PASI (outcome 11) rides v_pasi as an arm-level column; use it as
  # the fallback when an arm has no week-0 abs_pasi row.
  bp_key <- paste(df$ref_id, df$arm_no, sep = "|")
  bp_fallback <- list(
    mean = setNames(df$baseline_pasi_mean, bp_key),
    sd   = setNames(df$baseline_pasi_sd,   bp_key)
  )
  b <- baseline_lookup(df, "abs_pasi_mean", "abs_pasi_sd",
                       fallback = bp_fallback)
  d <- derive_change(b$mean, df$abs_pasi_mean, df$abs_pasi_change_mean,
                     baseline_sd = b$sd, follow_sd = df$abs_pasi_sd)
  df$baseline        <- fmt_mean_sd(b$mean, b$sd)
  df$on_tx           <- fmt_mean_sd(df$abs_pasi_mean, df$abs_pasi_sd)
  df$abs_pasi_change <- fmt_mean_sd_marked(d$change, df$abs_pasi_change_sd,
                                           d$change_derived,
                                           derived_sd = d$change_sd_derived)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  # Drop pure baseline rows; each follow-up row now carries its baseline.
  df <- df[is.na(df$timepoint) | df$timepoint > 0, ]
  has_any <- nzchar(df$on_tx) | nzchar(df$abs_pasi_change)
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", "baseline", "on_tx", "abs_pasi_change")]
}

format_dlqi_zero <- function(df) {
  format_binary_subset(df, c("dlqi_0_1", "dlqi_0"))
}

format_dlqi_absolute <- function(df) {
  b <- baseline_lookup(df, "abs_dlqi_mean", "abs_dlqi_sd")
  d <- derive_change(b$mean, df$abs_dlqi_mean, df$abs_dlqi_change_mean,
                     baseline_sd = b$sd, follow_sd = df$abs_dlqi_sd)
  df$baseline        <- fmt_mean_sd(b$mean, b$sd)
  df$on_tx           <- fmt_mean_sd(df$abs_dlqi_mean, df$abs_dlqi_sd)
  df$abs_dlqi_change <- fmt_mean_sd_marked(d$change, df$abs_dlqi_change_sd,
                                           d$change_derived,
                                           derived_sd = d$change_sd_derived)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  df <- df[is.na(df$timepoint) | df$timepoint > 0, ]
  has_any <- nzchar(df$on_tx) | nzchar(df$abs_dlqi_change)
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
        table    = "v_pasi",
        fmt      = format_pasi_absolute,
        colnames = c("Trial", "Drug", "N",
                     "Baseline", "Follow-up", "Δ from baseline"),
        note     = paste("Values marked with * are derived: the mean is",
                         "follow-up − baseline; the SD is approximated",
                         "assuming a baseline/follow-up correlation of 0.5.",
                         "Where an SD could not be derived it is omitted.")
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
      absolute = list(
        label    = "Absolute DLQI",
        table    = "v_dlqi",
        fmt      = format_dlqi_absolute,
        colnames = c("Trial", "Drug", "N",
                     "Baseline", "Follow-up", "Δ from baseline"),
        note     = paste("Values marked with * are derived: the mean is",
                         "follow-up − baseline; the SD is approximated",
                         "assuming a baseline/follow-up correlation of 0.5.",
                         "Where an SD could not be derived it is omitted.")
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
    zero     = function(df) survivors_binary(df,  c("dlqi_0_1", "dlqi_0")),
    absolute = function(df) survivors_absolute(df,
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
      list(code = "pasi50",  label = "PASI 50"),
      list(code = "pasi75",  label = "PASI 75"),
      list(code = "pasi90",  label = "PASI 90"),
      list(code = "pasi100", label = "PASI 100")
    ),
    absolute = list(
      list(code = "abs_pasi_change", label = "Δ from baseline (absolute PASI)")
    )
  ),
  dlqi = list(
    zero = list(
      list(code = "dlqi_0_1", label = "DLQI 0/1"),
      list(code = "dlqi_0",   label = "DLQI 0")
    ),
    absolute = list(
      list(code = "abs_dlqi_change", label = "Δ from baseline (absolute DLQI)")
    )
  ),
  safety = list(
    sae                = list(list(code = "sae", label = "Any SAE")),
    disc               = list(
      list(code = "disc_any", label = "Discontinuation (any)"),
      list(code = "disc_ae",  label = "Discontinuation (AE)")
    ),
    serious_infection  = list(
      list(code = "serious_infection", label = "Serious infection")
    ),
    injection_site_rxn = list(
      list(code = "injection_site_rxn", label = "Injection-site reaction")
    ),
    malignancy = list(
      list(code = "malignancy", label = "Malignancy")
    )
  )
)

# Single source of truth for whether the MA build step has been run.
ma_tables_present <- function() {
  con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  all(c("meta_analysis", "trial_estimates") %in% dbListTables(con))
}
HAS_MA <- ma_tables_present()

# Translate one MA "outcome spec" + the current filter into a rows/pooled
# bundle for forest_svg. Returns NULL when no data is available.
build_forest_inputs <- function(state, tab_id, outcome, response_method = "binomial") {
  is_continuous <- grepl("^abs_", outcome$code)
  is_harm       <- identical(tab_id, "safety")
  is_response   <- !is_continuous && !is_harm && tab_id %in% c("pasi", "dlqi")

  if (is.null(state)) {
    # No filter: one row per drug from the network model.
    # Comparisons vs Placebo may be stored in either orientation; collect both
    # and flip signs where Placebo is comp_tx. Falls back to the other effects
    # model when the preferred one has no data (e.g. rare-event endpoints with
    # only a fixed-effects network).
    measure_val <- if (is_continuous) "diff_cfb" else "rd"
    fetch_network_vs_placebo <- function(eff, method = NULL) {
      r_direct  <- fetch_ma(outcome$code, type = "network", effects = eff,
                            measure = measure_val, ref_tx = "Placebo", method = method)
      r_direct  <- r_direct[r_direct$comp_tx != "Placebo", , drop = FALSE]
      r_flipped <- fetch_ma(outcome$code, type = "network", effects = eff,
                            measure = measure_val, comp_tx = "Placebo", method = method)
      r_flipped <- r_flipped[r_flipped$ref_tx != "Placebo" &
                              !is.na(r_flipped$ref_tx) & r_flipped$ref_tx != "",
                             , drop = FALSE]
      if (nrow(r_flipped)) {
        tmp               <- r_flipped$lower
        r_flipped$lower   <- -r_flipped$upper
        r_flipped$upper   <- -tmp
        r_flipped$mean    <- -r_flipped$mean
        r_flipped$comp_tx <- r_flipped$ref_tx
        r_flipped$ref_tx  <- "Placebo"
      }
      rbind(r_direct, r_flipped)
    }

    if (is_response) {
      # Response endpoints: show FE + RE for the selected method (toggled in UI).
      method_lbl <- if (identical(response_method, "multinomial")) "multinomial" else "binomial"
      df_re <- fetch_network_vs_placebo("random", method_lbl)
      df_fe <- fetch_network_vs_placebo("fixed",  method_lbl)

      if (!nrow(df_re)) {
        return(list(empty_reason =
          sprintf("No %s random effects network model available for this endpoint.", method_lbl)))
      }

      common <- intersect(df_re$comp_tx, df_fe$comp_tx)
      df_re  <- df_re[match(common, df_re$comp_tx), , drop = FALSE]
      df_fe  <- df_fe[match(common, df_fe$comp_tx), , drop = FALSE]
      ord    <- order(df_re$mean, decreasing = TRUE)
      df_re  <- df_re[ord, , drop = FALSE]
      df_fe  <- df_fe[ord, , drop = FALSE]
      n_drugs <- nrow(df_re)

      PAIR_H   <- 14L
      DRUG_GAP <- 8L
      fe_gaps  <- c(0L, rep(DRUG_GAP, n_drugs - 1L))
      fe_tt <- mapply(function(drug, est, lo, hi)
        ma_tooltip(sprintf("%s — %s network estimate, fixed effects", drug, method_lbl), est, lo, hi, digits = 2),
        df_fe$comp_tx, df_fe$mean, df_fe$lower, df_fe$upper)
      re_tt <- mapply(function(drug, est, lo, hi)
        ma_tooltip(sprintf("%s — %s network estimate, random effects", drug, method_lbl), est, lo, hi, digits = 2),
        df_re$comp_tx, df_re$mean, df_re$lower, df_re$upper)
      rows <- data.frame(
        label     = c(rbind(df_fe$comp_tx, rep("", n_drugs))),
        badge     = c(rbind(rep("FE", n_drugs), rep("RE", n_drugs))),
        est       = c(rbind(df_fe$mean,  df_re$mean)),
        lo        = c(rbind(df_fe$lower, df_re$lower)),
        hi        = c(rbind(df_fe$upper, df_re$upper)),
        square_n  = NA_real_,
        klass     = c(rbind(rep("ma-square-fe", n_drugs), rep("ma-square-re", n_drugs))),
        row_h     = PAIR_H,
        gap_above = c(rbind(fe_gaps, rep(0L, n_drugs))),
        tooltip   = c(rbind(fe_tt, re_tt)),
        stringsAsFactors = FALSE
      )
    } else {
      df_re <- fetch_network_vs_placebo("random")
      if (!nrow(df_re)) {
        return(list(empty_reason =
          "No random effects network model available for this endpoint."))
      }
      df_fe <- fetch_network_vs_placebo("fixed")

      # Keep only drugs present in both models, sort by RE estimate best-first.
      common <- intersect(df_re$comp_tx, df_fe$comp_tx)
      df_re  <- df_re[match(common, df_re$comp_tx), , drop = FALSE]
      df_fe  <- df_fe[match(common, df_fe$comp_tx), , drop = FALSE]
      ord    <- order(df_re$mean, decreasing = !(is_continuous || is_harm))
      df_re  <- df_re[ord, , drop = FALSE]
      df_fe  <- df_fe[ord, , drop = FALSE]
      n_drugs <- nrow(df_re)

      # Interleave FE (grey) then RE (blue) per drug; add gap between drug pairs.
      PAIR_H   <- 14L
      DRUG_GAP <- 8L
      fe_gaps  <- c(0L, rep(DRUG_GAP, n_drugs - 1L))
      fe_tt <- mapply(function(drug, est, lo, hi)
        ma_tooltip(sprintf("%s — network estimate, fixed effects", drug), est, lo, hi, digits = 2),
        df_fe$comp_tx, df_fe$mean, df_fe$lower, df_fe$upper)
      re_tt <- mapply(function(drug, est, lo, hi)
        ma_tooltip(sprintf("%s — network estimate, random effects", drug), est, lo, hi, digits = 2),
        df_re$comp_tx, df_re$mean, df_re$lower, df_re$upper)
      rows <- data.frame(
        label     = c(rbind(df_fe$comp_tx, rep("", n_drugs))),
        badge     = c(rbind(rep("FE", n_drugs), rep("RE", n_drugs))),
        est       = c(rbind(df_fe$mean,  df_re$mean)),
        lo        = c(rbind(df_fe$lower, df_re$lower)),
        hi        = c(rbind(df_fe$upper, df_re$upper)),
        square_n  = NA_real_,
        klass     = c(rbind(rep("ma-square-fe", n_drugs),
                            rep("ma-square-re", n_drugs))),
        row_h     = PAIR_H,
        gap_above = c(rbind(fe_gaps, rep(0L, n_drugs))),
        tooltip   = c(rbind(fe_tt, re_tt)),
        stringsAsFactors = FALSE
      )
    }

    axis_label <- if (is_continuous)
      "Difference in change from baseline vs Placebo (95% CI)"
    else
      "Risk difference vs Placebo (95% CI)"

    if (is_continuous || is_harm) {
      dir_left  <- "← favours drug"
      dir_right <- "favours placebo →"
    } else {
      dir_left  <- "← favours placebo"
      dir_right <- "favours drug →"
    }

    return(list(rows = rows, pooled = NULL,
                effective_scale = "md",
                axis_label = axis_label,
                dir_left = dir_left, dir_right = dir_right,
                n_drugs = n_drugs,
                response_method = if (is_response) method_lbl else NULL))
  }

  if (identical(state$kind, "edge")) {
    # Edge filter: per-trial pairwise forest with FE/RE + network diamonds.
    measure_val <- if (is_continuous) "diff_cfb" else "rd"

    # Determine canonical orientation: pairwise MA preferred (it may not exist
    # for single-trial drugs), trial estimates as fallback. Both (A,B) and
    # (B,A) are tried in each case.
    find_orient <- function(a, b) {
      r <- fetch_ma(outcome$code, type = "pairwise", effects = "fixed",
                    comp_tx = a, ref_tx = b, measure = measure_val)
      if (nrow(r)) return(c(a, b))
      r <- fetch_trials(outcome$code, comp_tx = a, ref_tx = b,
                        measure = measure_val)
      if (nrow(r)) return(c(a, b))
      if (is_response) {
        r <- fetch_ma(outcome$code, type = "network", effects = "random",
                      comp_tx = a, ref_tx = b, measure = measure_val,
                      method = "multinomial")
        if (nrow(r)) return(c(a, b))
      }
      NULL
    }
    orient <- find_orient(state$from, state$to)
    if (is.null(orient)) orient <- find_orient(state$to, state$from)
    if (is.null(orient)) return(NULL)
    comp <- orient[1]; ref <- orient[2]

    trials <- fetch_trials(outcome$code, comp_tx = comp, ref_tx = ref,
                           measure = measure_val)
    if (!nrow(trials) && !is_response) return(NULL)

    if (nrow(trials)) {
      # When one side is Placebo, normalise to drug-vs-Placebo (positive = favours drug).
      if (identical(comp, "Placebo") && !identical(ref, "Placebo")) {
        trials$mean  <- -trials$mean
        tmp          <- trials$lower
        trials$lower <- -trials$upper
        trials$upper <- -tmp
        for (pair in list(c("comp_tx", "ref_tx"), c("n_tx", "n_ref"), c("k_tx", "k_ref"),
                          c("mean_tx", "mean_ref"), c("sd_tx", "sd_ref"))) {
          tmp               <- trials[[pair[1]]]
          trials[[pair[1]]] <- trials[[pair[2]]]
          trials[[pair[2]]] <- tmp
        }
        comp <- trials$comp_tx[1]
        ref  <- trials$ref_tx[1]
      }

      rows <- data.frame(
        label    = trial_name[as.character(trials$ref_id)],
        est      = trials$mean, lo = trials$lower, hi = trials$upper,
        square_n = trials$n_tx + coalesce0(trials$n_ref),
        stringsAsFactors = FALSE
      )
      rows$tooltip <- vapply(seq_len(nrow(rows)), function(i) {
        t <- trials[i, ]
        lbl_tx  <- sprintf("%s (events/n)", t$comp_tx)
        lbl_ref <- sprintf("%s (events/n)", t$ref_tx)
        if (!is_continuous)
          ma_tooltip(trial_name[as.character(t$ref_id)] %||% as.character(t$ref_id), t$mean, t$lower, t$upper,
                     extra = setNames(
                       c(sprintf("%d / %d", t$k_tx, t$n_tx),
                         sprintf("%d / %d", t$k_ref, t$n_ref)),
                       c(lbl_tx, lbl_ref)))
        else
          ma_tooltip(trial_name[as.character(t$ref_id)] %||% as.character(t$ref_id), t$mean, t$lower, t$upper,
                     extra = setNames(
                       c(sprintf("%.2f (%.2f), %d", t$mean_tx, t$sd_tx, t$n_tx),
                         sprintf("%.2f (%.2f), %d", t$mean_ref, t$sd_ref, t$n_ref)),
                       c(sprintf("%s: mean (SD), n", t$comp_tx),
                         sprintf("%s: mean (SD), n", t$ref_tx))))
      }, character(1))
    } else {
      if (identical(comp, "Placebo") && !identical(ref, "Placebo")) {
        tmp  <- comp; comp <- ref; ref <- tmp
      }
      rows <- data.frame(label = character(), est = numeric(), lo = numeric(),
                         hi = numeric(), square_n = numeric(), tooltip = character(),
                         stringsAsFactors = FALSE)
    }

    # Pooled diamonds: pairwise FE/RE + network FE/RE (all in comp→ref direction).
    pooled_list <- list()
    for (kind_pair in list(c("FE", "fixed"), c("RE", "random"))) {
      r <- fetch_ma_directed(outcome$code, "pairwise", kind_pair[2],
                             comp, ref, measure_val)
      if (!is.null(r))
        pooled_list[[kind_pair[1]]] <- data.frame(
          kind = kind_pair[1], est = r$mean, lo = r$lower, hi = r$upper)
    }
    if (is_response) {
      for (kind_pair in list(c("Bin-NMA-FE",  "fixed",  "binomial"),
                             c("Bin-NMA-RE",  "random", "binomial"),
                             c("Mult-NMA-FE", "fixed",  "multinomial"),
                             c("Mult-NMA-RE", "random", "multinomial"))) {
        r <- fetch_ma_directed(outcome$code, "network", kind_pair[2],
                               comp, ref, measure_val, method = kind_pair[3])
        if (!is.null(r))
          pooled_list[[kind_pair[1]]] <- data.frame(
            kind = kind_pair[1], est = r$mean, lo = r$lower, hi = r$upper)
      }
    } else {
      for (kind_pair in list(c("NMA-FE", "fixed"), c("NMA-RE", "random"))) {
        r <- fetch_ma_directed(outcome$code, "network", kind_pair[2],
                               comp, ref, measure_val)
        if (!is.null(r))
          pooled_list[[kind_pair[1]]] <- data.frame(
            kind = kind_pair[1], est = r$mean, lo = r$lower, hi = r$upper)
      }
    }
    pooled <- if (length(pooled_list)) do.call(rbind, pooled_list) else NULL
    if (!nrow(rows) && is.null(pooled)) return(NULL)

    dir_left <- dir_right <- NULL
    if (!is_continuous && is_harm) {
      dir_left  <- sprintf("← favours %s", comp)
      dir_right <- sprintf("favours %s →", ref)
    } else if (!is_continuous) {
      dir_left  <- sprintf("← favours %s", ref)
      dir_right <- sprintf("favours %s →", comp)
    } else {
      # Continuous: lower CFB difference favours comp (less change = better?
      # Actually lower PASI/DLQI is better so more negative diff_cfb favours comp).
      dir_left  <- sprintf("← favours %s", comp)
      dir_right <- sprintf("favours %s →", ref)
    }

    return(list(rows = rows, pooled = pooled,
                effective_scale = "md",
                comparison = sprintf("%s vs %s", comp, ref),
                axis_label = if (!is_continuous)
                  "Risk difference (95% CI)"
                else
                  "Mean difference in change from baseline (95% CI)",
                dir_left = dir_left, dir_right = dir_right))
  }

  if (identical(state$kind, "node")) {
    # Node filter: per-trial single-arm estimates with univariate pooled.
    measure_val <- if (is_continuous) "cfb" else "rate"
    eff_scale   <- if (is_continuous) "md" else "prop"

    trials <- fetch_trials(outcome$code, comp_tx = state$drug,
                           measure = measure_val)
    if (!nrow(trials) && !is_response) return(NULL)

    digits <- if (eff_scale == "prop") 3 else 2
    if (nrow(trials)) {
      rows <- data.frame(
        label    = trial_name[as.character(trials$ref_id)],
        est      = trials$mean, lo = trials$lower, hi = trials$upper,
        square_n = trials$n_tx,
        stringsAsFactors = FALSE
      )
      event_lbl <- if (is_harm) "Events" else "Responders"
      rows$tooltip <- vapply(seq_len(nrow(rows)), function(i) {
        t <- trials[i, ]
        if (!is_continuous)
          ma_tooltip(trial_name[as.character(t$ref_id)] %||% as.character(t$ref_id), t$mean, t$lower, t$upper,
                     extra = c(setNames(sprintf("%d / %d", t$k_tx, t$n_tx),
                                        event_lbl)),
                     digits = digits)
        else
          ma_tooltip(trial_name[as.character(t$ref_id)] %||% as.character(t$ref_id), t$mean, t$lower, t$upper,
                     extra = c("Mean (SD), n" = sprintf("%.2f (%.2f), %d",
                                                         t$mean, t$sd_tx, t$n_tx)),
                     digits = digits)
      }, character(1))
    } else {
      rows <- data.frame(label = character(), est = numeric(), lo = numeric(),
                         hi = numeric(), square_n = numeric(), tooltip = character(),
                         stringsAsFactors = FALSE)
    }

    pooled_list <- list()
    for (kind_pair in list(c("Pool-FE", "fixed"), c("Pool-RE", "random"))) {
      r <- fetch_ma(outcome$code, type = "univariate", effects = kind_pair[2],
                    comp_tx = state$drug, measure = measure_val)
      if (nrow(r))
        pooled_list[[kind_pair[1]]] <- data.frame(
          kind = kind_pair[1], est = r$mean[1], lo = r$lower[1], hi = r$upper[1])
    }
    if (is_response) {
      for (kind_pair in list(c("Bin-NMA-R-FE",  "fixed",  "binomial"),
                             c("Bin-NMA-R-RE",  "random", "binomial"),
                             c("Mult-NMA-R-FE", "fixed",  "multinomial"),
                             c("Mult-NMA-R-RE", "random", "multinomial"))) {
        r <- fetch_ma(outcome$code, type = "network", effects = kind_pair[2],
                      comp_tx = state$drug, measure = measure_val, method = kind_pair[3])
        if (nrow(r))
          pooled_list[[kind_pair[1]]] <- data.frame(
            kind = kind_pair[1], est = r$mean[1], lo = r$lower[1], hi = r$upper[1])
      }
    } else {
      for (kind_pair in list(c("NMA-R-FE", "fixed"), c("NMA-R-RE", "random"))) {
        r <- fetch_ma(outcome$code, type = "network", effects = kind_pair[2],
                      comp_tx = state$drug, measure = measure_val)
        if (nrow(r))
          pooled_list[[kind_pair[1]]] <- data.frame(
            kind = kind_pair[1], est = r$mean[1], lo = r$lower[1], hi = r$upper[1])
      }
    }
    pooled <- if (length(pooled_list)) do.call(rbind, pooled_list) else NULL
    if (!nrow(rows) && is.null(pooled)) return(NULL)

    prop_lbl <- if (is_harm) "Event rate" else "Response rate"
    axis_label <- if (is_continuous)
      sprintf("Change from baseline (95%% CI), %s", state$drug)
    else
      sprintf("%s, %s (95%% CI)", prop_lbl, state$drug)

    return(list(rows = rows, pooled = pooled,
                effective_scale = eff_scale,
                drug = state$drug,
                axis_label = axis_label,
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

    // Trial-name click handler. Sends the ref-id to Shiny so the server
    // can open a detail modal. Event delegation survives DT redraws.
    document.addEventListener('click', function(ev) {
      var t = ev.target.closest && ev.target.closest('.trial-pop');
      if (t) {
        ev.preventDefault();
        Shiny.setInputValue('trial_modal', t.getAttribute('data-ref-id'),
                            {priority: 'event'});
      }
    });

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
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com",
              crossorigin = NA),
    tags$link(rel = "stylesheet",
              href = paste0("https://fonts.googleapis.com/css2?",
                            "family=Inter:wght@400;500;600;700&display=swap"))
  ),
  tags$head(tags$link(rel = "stylesheet", href = "style.css")),
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
        click = "function(p){
          if(p.nodes && p.nodes.length)      Shiny.setInputValue('nma_node',  p.nodes[0],      {priority:'event'});
          else if(p.edges && p.edges.length) Shiny.setInputValue('nma_edge',  p.edges[0],      {priority:'event'});
          else                               Shiny.setInputValue('nma_clear', Math.random(),   {priority:'event'});
        }"
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
                           group_lbl = NULL, is_response_nofilter = FALSE,
                           response_method = "multinomial")

  # Build a single plot block for one outcome under the active state.
  build_plot_block <- function(state, tab_id, outc, response_method = "binomial") {
    inputs <- tryCatch(build_forest_inputs(state, tab_id, outc, response_method),
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
    eff_scale <- inputs$effective_scale %||% "md"
    stats <- if (!is.null(inputs$pooled) && nrow(inputs$pooled) > 0) {
      ax_lbl <- inputs$axis_label %||% switch(eff_scale, md = "MD", prop = "Proportion", "Effect")
      n_st <- nrow(inputs$rows)
      sprintf("%s • %d %s", ax_lbl, n_st, if (n_st == 1) "study" else "studies")
    } else if (!is.null(inputs$n_drugs)) {
      if (!is.null(inputs$response_method))
        sprintf("Network meta-analysis, FE + RE — %d drugs (%s)", inputs$n_drugs, inputs$response_method)
      else
        sprintf("Network meta-analysis, FE + RE — %d drugs", inputs$n_drugs)
    } else {
      sprintf("%d %s", nrow(inputs$rows),
              if (is.null(state)) "drugs" else "studies")
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

  # Reactively rendered plot area. Driven by ma_ctx (set on modal open).
  output$ma_plot_area <- renderUI({
    req(ma_ctx$active, ma_ctx$tab_id, ma_ctx$gid, ma_ctx$outcomes)
    tab_id          <- ma_ctx$tab_id
    state           <- ma_ctx$state
    outcomes        <- ma_ctx$outcomes
    response_method <- ma_ctx$response_method
    lapply(outcomes, function(outc) build_plot_block(state, tab_id, outc, response_method))
  })

  observeEvent(input$ma_method, {
    ma_ctx$response_method <- input$ma_method
  })

  output$ma_method_toggle <- renderUI({
    req(ma_ctx$is_response_nofilter)
    div(class = "ma-method-toggle",
        radioButtons("ma_method", label = NULL,
                     choices = c("Binomial" = "binomial", "Multinomial" = "multinomial"),
                     selected = "multinomial",
                     inline = TRUE))
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
    ma_ctx$tab_id               <- tab_id
    ma_ctx$state                <- state
    ma_ctx$gid                  <- gid
    ma_ctx$outcomes             <- outcomes
    ma_ctx$state_lbl            <- state_lbl
    ma_ctx$group_lbl            <- group_lbl
    ma_ctx$is_response_nofilter <- is.null(state) && tab_id %in% c("pasi", "dlqi") &&
                                   gid %in% c("response", "zero")
    ma_ctx$response_method      <- "multinomial"
    ma_ctx$active               <- TRUE

    summary_text <- sprintf("%s | endpoint group: %s", state_lbl, group_lbl)

    showModal(modalDialog(
      title = tagList(icon("chart-column"), " Meta-analysis"),
      easyClose = TRUE, size = "l", class = "ma-modal",
      footer = modalButton("Close"),
      div(
        div(class = "ma-summary", summary_text),
        uiOutput("ma_method_toggle"),
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

  # --- Trial detail modal ---------------------------------------------------
  observeEvent(input$trial_modal, {
    sid  <- input$trial_modal
    name <- trial_name[sid] %||% sid

    # Primary publication title, shown as a subtitle under the trial name.
    pubs <- pubs_by_study[[sid]]
    primary_title <- if (!is.null(pubs) && any(pubs$is_primary == 1)) {
      trimws(pubs$title[which(pubs$is_primary == 1)[1]])
    } else NA_character_
    modal_title <- if (!is.na(primary_title) && nzchar(primary_title)) {
      tagList(
        div(class = "trial-modal-title", name),
        div(class = "trial-modal-subtitle", primary_title)
      )
    } else {
      name
    }

    # References
    refs_html <- cites_html_by_study[[sid]] %||%
      '<p class="trial-modal-empty">No references available.</p>'

    # Baseline characteristics
    bl <- fetch_baselines(as.integer(sid))
    if (nrow(bl)) {
      bl$subcategory[bl$subcategory == "Comorbidity"] <- "Psoriasis characteristics"

      arms <- unique(bl[, c("arm_no", "drug", "dose")])
      arms <- arms[order(arms$arm_no), ]
      # Header: drug on top (spanning same-drug arms), dose (n = ...) below.
      arm_n <- tapply(bl$n, bl$arm_no, function(x) max(x, na.rm = TRUE))
      has_dose <- !is.na(arms$dose) & nzchar(arms$dose)
      dose_label <- ifelse(has_dose, arms$dose, "")
      n_vals <- arm_n[as.character(arms$arm_no)]
      has_n <- !is.na(n_vals) & is.finite(n_vals)
      dose_label[has_n] <- trimws(sprintf("%s (n = %d)", dose_label[has_n],
                                          as.integer(n_vals[has_n])))

      bl_header <- build_arm_header(arms$drug, dose_label)

      subcat_order <- c("Demographics", "Psoriasis characteristics",
                        "Previous therapy")
      subcats <- intersect(subcat_order, unique(bl$subcategory))
      bl_rows <- character(0)
      for (sc in subcats) {
        sc_rows <- character(0)
        labels <- unique(bl$label[bl$subcategory == sc])
        for (lbl in labels) {
          cells <- character(nrow(arms))
          for (j in seq_len(nrow(arms))) {
            row <- bl[bl$subcategory == sc & bl$label == lbl &
                      bl$arm_no == arms$arm_no[j], ]
            if (!nrow(row)) { cells[j] <- ""; next }
            r <- row[1, ]
            if (!is.na(r$k) && !is.na(r$n) && r$n > 0) {
              cells[j] <- sprintf("%d (%d%%)", as.integer(r$k),
                                  round(r$k / r$n * 100))
            } else if (!is.na(r$mean)) {
              cells[j] <- fmt_mean_sd(r$mean, r$sd)
            } else {
              cells[j] <- ""
            }
          }
          if (all(!nzchar(cells))) next
          sc_rows <- c(sc_rows, paste0(
            "<tr><td>", htmltools::htmlEscape(lbl), "</td>",
            paste0("<td>", cells, "</td>", collapse = ""),
            "</tr>"))
        }
        if (!length(sc_rows)) next
        bl_rows <- c(bl_rows,
          sprintf('<tr class="subcat-row"><td colspan="%d">%s</td></tr>',
                  1L + nrow(arms), htmltools::htmlEscape(sc)),
          sc_rows)
      }
      bl_html <- paste0(
        "<h4>Baseline characteristics</h4>",
        "<table>", bl_header, paste(bl_rows, collapse = ""), "</table>")
    } else {
      bl_html <- '<p class="trial-modal-empty">No baseline data available.</p>'
    }

    # Results from all views
    res <- fetch_trial_results(as.integer(sid))
    results_html <- character(0)

    view_labels <- list(
      pasi = list(
        title = "PASI",
        binary = c(pasi50 = "PASI 50", pasi75 = "PASI 75",
                   pasi90 = "PASI 90", pasi100 = "PASI 100"),
        continuous = c(abs_pasi_mean = "Absolute PASI")
      ),
      dlqi = list(
        title = "DLQI",
        binary = c(dlqi_0_1 = "DLQI 0/1", dlqi_0 = "DLQI 0"),
        continuous = c(abs_dlqi_mean = "Absolute DLQI")
      ),
      safety = list(
        title = "Safety",
        binary = c(sae = "Any SAE", disc_any = "Disc. (any)",
                   disc_ae = "Disc. (AE)",
                   serious_infection = "Serious infection",
                   injection_site_rxn = "Injection site rxn",
                   malignancy = "Malignancy"),
        continuous = character(0)
      )
    )

    for (vname in names(view_labels)) {
      df <- res[[vname]]
      if (!nrow(df)) next
      vl <- view_labels[[vname]]

      # Arms that have data in this view — header: drug on top, dose (n = ...) below.
      arm_cols <- unique(df[, c("arm_no", "drug", "dose")])
      arm_cols <- arm_cols[order(arm_cols$arm_no), ]
      has_dose <- !is.na(arm_cols$dose) & nzchar(arm_cols$dose)
      dose_label <- ifelse(has_dose, arm_cols$dose, "")
      res_arm_n <- tapply(df$n, df$arm_no, function(x) max(x, na.rm = TRUE))
      n_vals <- res_arm_n[as.character(arm_cols$arm_no)]
      has_n <- !is.na(n_vals) & is.finite(n_vals)
      dose_label[has_n] <- trimws(sprintf("%s (n = %d)", dose_label[has_n],
                                          as.integer(n_vals[has_n])))

      header <- build_arm_header(arm_cols$drug, dose_label)

      timepoints <- sort(unique(df$timepoint[!is.na(df$timepoint) &
                                             df$timepoint > 0]))
      if (!length(timepoints)) timepoints <- unique(df$timepoint)

      rows <- character(0)
      # Binary outcomes
      for (col in names(vl$binary)) {
        if (!(col %in% names(df))) next
        if (all(is.na(df[[col]]))) next
        for (tp in timepoints) {
          tp_df <- df[!is.na(df$timepoint) & df$timepoint == tp, ]
          if (!nrow(tp_df)) next
          cells <- character(nrow(arm_cols))
          for (j in seq_len(nrow(arm_cols))) {
            r <- tp_df[tp_df$arm_no == arm_cols$arm_no[j], ]
            if (!nrow(r) || is.na(r[[col]][1])) { cells[j] <- ""; next }
            cells[j] <- fmt_pasi(r[[col]][1], r$n[1])
          }
          if (all(!nzchar(cells))) next
          tp_unit <- df$timepoint_unit[1]
          if (!is.na(tp_unit) && !grepl("s$", tp_unit)) tp_unit <- paste0(tp_unit, "s")
          tp_label <- paste0(tp, " ", tp_unit)
          row_label <- sprintf("%s (%s)", vl$binary[[col]], tp_label)
          rows <- c(rows, paste0(
            "<tr><td>", htmltools::htmlEscape(row_label), "</td>",
            paste0("<td>", cells, "</td>", collapse = ""),
            "</tr>"))
        }
      }
      # Continuous outcomes
      for (col in names(vl$continuous)) {
        sd_col <- sub("_mean$", "_sd", col)
        if (!(col %in% names(df))) next
        if (all(is.na(df[[col]]))) next
        for (tp in timepoints) {
          tp_df <- df[!is.na(df$timepoint) & df$timepoint == tp, ]
          if (!nrow(tp_df)) next
          cells <- character(nrow(arm_cols))
          for (j in seq_len(nrow(arm_cols))) {
            r <- tp_df[tp_df$arm_no == arm_cols$arm_no[j], ]
            if (!nrow(r) || is.na(r[[col]][1])) { cells[j] <- ""; next }
            sd_val <- if (sd_col %in% names(r)) r[[sd_col]][1] else NA
            cells[j] <- fmt_mean_sd(r[[col]][1], sd_val)
          }
          if (all(!nzchar(cells))) next
          tp_unit <- df$timepoint_unit[1]
          if (!is.na(tp_unit) && !grepl("s$", tp_unit)) tp_unit <- paste0(tp_unit, "s")
          tp_label <- paste0(tp, " ", tp_unit)
          row_label <- sprintf("%s (%s)", vl$continuous[[col]], tp_label)
          rows <- c(rows, paste0(
            "<tr><td>", htmltools::htmlEscape(row_label), "</td>",
            paste0("<td>", cells, "</td>", collapse = ""),
            "</tr>"))
        }
      }

      if (length(rows)) {
        results_html <- c(results_html, paste0(
          "<h4>", htmltools::htmlEscape(vl$title), "</h4>",
          "<table>", header, paste(rows, collapse = ""), "</table>"))
      }
    }

    if (!length(results_html)) {
      results_html <- '<p class="trial-modal-empty">No results available.</p>'
    }

    showModal(modalDialog(
      title = modal_title,
      size  = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      div(class = "trial-modal-layout",
        div(class = "trial-modal-refs", HTML(refs_html)),
        div(class = "trial-modal-data",
            HTML(bl_html),
            HTML(paste(results_html, collapse = "")))
      )
    ) |> tagAppendAttributes(class = "trial-modal"))
  })

  output$download_db <- downloadHandler(
    filename    = function() "psoriasis-rcts.sqlite",
    contentType = "application/x-sqlite3",
    content     = function(file) file.copy(DB_PATH, file)
  )
}

shinyApp(ui, server)
