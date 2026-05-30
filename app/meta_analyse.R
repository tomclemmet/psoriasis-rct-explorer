# Offline meta-analysis DRIVER. Reads psoriasis-rcts.sqlite, prepares the data
# for every pooled estimate the Shiny app might display, hands each prepared
# data frame to a model function in ma_models.R, and writes the returned results
# back into the same SQLite file in six tables:
#
#   ma_pairwise          - one row per (drug_a, drug_b, endpoint, outcome)
#                          comparison (FE + RE pooled estimate + heterogeneity).
#   ma_pairwise_trials   - per-trial contributions to each ma_pairwise row.
#   ma_proportion        - one row per (drug, endpoint, outcome) pooled response
#                          rate (single-arm meta-analysis).
#   ma_proportion_trials - per-arm contributions to each ma_proportion row.
#   ma_nma               - one network-level summary per (endpoint, outcome).
#   ma_nma_estimates     - one row per ordered drug pair from each network.
#
# DIVISION OF LABOUR
#   This file (the driver / plumbing) owns: reading the views, primary-timepoint
#   selection, change derivation, dose aggregation, comparison enumeration, and
#   the SQLite write. It does NOT fit any models.
#   ma_models.R (the file you edit) owns the statistics: it receives a prepared
#   data frame per analysis unit and returns results in the fixed contract. Swap
#   in JAGS/STAN/brms there without touching this file or app.R.
#
# Per-trial primary timepoint: the latest measurement at or before a target week
# per (trial, arm, outcome). Targets: 16 wk for PASI, 24 wk for DLQI and safety.
#
# Outcomes are grouped into "specs" (see the `analyses` registry below). The
# default models fit each outcome in a spec independently -- identical to the
# original per-outcome analysis. A joint model (e.g. ordinal PASI) lives in one
# spec whose model fits all its outcomes at once; this file needs no change for
# that, because every result row carries its own outcome_code.
#
# USAGE
#   Rscript app/meta_analyse.R                 # rebuild all families
#   Rscript app/meta_analyse.R --only nma      # rebuild only ma_nma*, leaving
#                                              # the other tables untouched
#   families: pairwise | proportion | nma  (comma-separated, e.g. --only nma,pairwise)
#
# Re-running overwrites the tables for the families that ran. Cheap (~hundreds
# of meta::* calls) for the default frequentist models.

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(meta)
  library(netmeta)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

resolve_script_dir <- function() {
  file_arg <- sub("^--file=", "",
                  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg)) return(normalizePath(dirname(file_arg[1])))
  of <- try(sys.frame(1)$ofile, silent = TRUE)
  if (!inherits(of, "try-error") && length(of) && nzchar(of))
    return(normalizePath(dirname(of)))
  normalizePath(".")
}
SCRIPT_DIR <- resolve_script_dir()
DB_PATH    <- file.path(SCRIPT_DIR, "psoriasis-rcts.sqlite")
if (!file.exists(DB_PATH)) stop("psoriasis-rcts.sqlite not found - run convert.R first.")

# The model functions + result contract live next door.
source(file.path(SCRIPT_DIR, "ma_models.R"))

# meta:::settings.meta() chatter is just noise here.
suppressMessages(settings.meta(CIbracket = "(", print.tau2 = FALSE))

# --only <family[,family]> selective rebuild ---------------------------------
ALL_FAMILIES <- c("pairwise", "proportion", "nma")
parse_only <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  val <- NULL
  i <- which(a == "--only")
  if (length(i) && i[1] < length(a)) val <- a[i[1] + 1L]
  eq <- grep("^--only=", a, value = TRUE)
  if (length(eq)) val <- sub("^--only=", "", eq[1])
  if (is.null(val)) return(ALL_FAMILIES)
  fams <- trimws(strsplit(val, ",")[[1]])
  bad <- setdiff(fams, ALL_FAMILIES)
  if (length(bad)) stop(sprintf("--only: unknown family '%s' (valid: %s)",
                                paste(bad, collapse = ", "),
                                paste(ALL_FAMILIES, collapse = ", ")))
  fams
}
RUN_FAMILIES <- parse_only()

# ---------------------------------------------------------------------------
# Analysis registry. Each spec groups outcomes that share a view, effect kind
# and primary-timepoint target. The default per-outcome models make grouping
# cosmetic (output is identical however outcomes are grouped); grouping exists
# so a future joint model can claim a whole set (e.g. all four PASI thresholds).
# Optional per-family overrides: model_pairwise / model_proportion / model_nma.
# ---------------------------------------------------------------------------
analyses <- list(
  list(group = "pasi", view = "v_pasi", kind = "binary", target_wk = 16,
       outcomes = list(
         list(code = "pasi50",  k_col = "pasi50"),
         list(code = "pasi75",  k_col = "pasi75"),
         list(code = "pasi90",  k_col = "pasi90"),
         list(code = "pasi100", k_col = "pasi100"))),
  list(group = "pasi_abs", view = "v_pasi_abs", kind = "continuous", target_wk = 16,
       outcomes = list(
         list(code = "abs_pasi_change",
              mean_col = "abs_pasi_change_mean", sd_col = "abs_pasi_change_sd",
              baseline_mean_col = "abs_pasi_mean", baseline_sd_col = "abs_pasi_sd",
              followup_mean_col = "abs_pasi_mean", followup_sd_col = "abs_pasi_sd"))),
  list(group = "dlqi", view = "v_dlqi", kind = "binary", target_wk = 24,
       outcomes = list(
         list(code = "dlqi_0_1",     k_col = "dlqi_0_1"),
         list(code = "dlqi_0",       k_col = "dlqi_0"),
         list(code = "dlqi_le5",     k_col = "dlqi_le5"),
         list(code = "dlqi_5pt_dec", k_col = "dlqi_5pt_dec"),
         list(code = "dlqi_4pt_dec", k_col = "dlqi_4pt_dec"))),
  list(group = "dlqi_abs", view = "v_dlqi", kind = "continuous", target_wk = 24,
       outcomes = list(
         list(code = "abs_dlqi_change",
              mean_col = "abs_dlqi_change_mean", sd_col = "abs_dlqi_change_sd",
              baseline_mean_col = "abs_dlqi_mean", baseline_sd_col = "abs_dlqi_sd",
              followup_mean_col = "abs_dlqi_mean", followup_sd_col = "abs_dlqi_sd"))),
  list(group = "safety", view = "v_safety", kind = "binary", target_wk = 24,
       outcomes = list(
         list(code = "sae",                k_col = "sae"),
         list(code = "disc_any",           k_col = "disc_any"),
         list(code = "disc_ae",            k_col = "disc_ae"),
         list(code = "serious_infection",  k_col = "serious_infection"),
         list(code = "injection_site_rxn", k_col = "injection_site_rxn"),
         list(code = "malignancy",         k_col = "malignancy"),
         list(code = "nmsc",               k_col = "nmsc"),
         list(code = "malignancy_non_nmsc", k_col = "malignancy_non_nmsc")))
)

# Per-outcome view of a spec, shaped like the old `endpoints` entries so the
# data-prep helpers below need no changes.
spec_ep <- function(spec, oc)
  c(list(group = spec$group, view = spec$view, kind = spec$kind,
         target_wk = spec$target_wk), oc)

# ---------------------------------------------------------------------------
# Pull each view once (a view may back several outcomes).
# ---------------------------------------------------------------------------
con <- dbConnect(SQLite(), DB_PATH)
on.exit(dbDisconnect(con), add = TRUE)

view_cache <- new.env(parent = emptyenv())
get_view <- function(name) {
  if (is.null(view_cache[[name]]))
    view_cache[[name]] <- dbGetQuery(con, sprintf("SELECT * FROM %s", name))
  view_cache[[name]]
}

# Per-arm primary timepoint: latest row at or before target_wk for that
# (ref_id, arm_no) with a non-NA value in `value_col`.
pick_primary_timepoint <- function(df, value_col, target_wk) {
  ok <- !is.na(df[[value_col]]) & !is.na(df$timepoint) & df$timepoint > 0 &
        df$timepoint <= target_wk
  d <- df[ok, , drop = FALSE]
  if (!nrow(d)) return(d[0, , drop = FALSE])
  key <- paste(d$ref_id, d$arm_no, sep = "|")
  ord <- order(key, -d$timepoint)
  d   <- d[ord, , drop = FALSE]
  d[!duplicated(paste(d$ref_id, d$arm_no, sep = "|")), , drop = FALSE]
}

# Continuous variant: an arm contributes if EITHER its direct change is reported
# OR baseline + follow-up are both present (we derive the change).
pick_primary_timepoint_continuous <- function(df, ep, target_wk) {
  has_change <- !is.na(df[[ep$mean_col]])
  has_bf     <- !is.na(df[[ep$followup_mean_col]])
  ok <- (has_change | has_bf) &
        !is.na(df$timepoint) & df$timepoint > 0 & df$timepoint <= target_wk
  d <- df[ok, , drop = FALSE]
  if (!nrow(d)) return(d[0, , drop = FALSE])
  key <- paste(d$ref_id, d$arm_no, sep = "|")
  ord <- order(key, -d$timepoint)
  d   <- d[ord, , drop = FALSE]
  d[!duplicated(paste(d$ref_id, d$arm_no, sep = "|")), , drop = FALSE]
}

# Fill change_mean = follow - baseline when only baseline + follow-up reported.
attach_derived_change <- function(d, raw, ep) {
  if (!nrow(d)) {
    d$change_mean <- numeric(0); d$change_sd <- numeric(0); return(d)
  }
  is_b <- !is.na(raw$timepoint) & raw$timepoint == 0
  bkey <- paste(raw$ref_id[is_b], raw$arm_no[is_b], sep = "|")
  bmean <- raw[[ep$baseline_mean_col]][is_b]
  bsd   <- raw[[ep$baseline_sd_col]][is_b]
  key <- paste(d$ref_id, d$arm_no, sep = "|")
  i <- match(key, bkey)
  baseline_mean <- bmean[i]
  baseline_sd   <- bsd[i]
  change_mean <- d[[ep$mean_col]]
  change_sd   <- d[[ep$sd_col]]
  derive <- is.na(change_mean) & !is.na(d[[ep$followup_mean_col]]) &
            !is.na(baseline_mean)
  change_mean[derive] <- d[[ep$followup_mean_col]][derive] - baseline_mean[derive]
  change_sd[derive]   <- NA_real_
  d$change_mean <- change_mean
  d$change_sd   <- change_sd
  d
}

# Per-endpoint "long" table of contributing arms: one row per
# (ref_id, trial, arm_no, drug, n, k_or_mean[, sd]).
build_arm_table <- function(ep) {
  raw <- get_view(ep$view)
  if (ep$kind == "binary") {
    d <- pick_primary_timepoint(raw, ep$k_col, ep$target_wk)
    if (!nrow(d)) return(data.frame())
    d <- d[!is.na(d$n) & d$n > 0 & !is.na(d$drug) & nzchar(d$drug), , drop = FALSE]
    data.frame(
      ref_id    = d$ref_id,
      trial     = d$trial,
      arm_no    = d$arm_no,
      drug      = d$drug,
      timepoint = d$timepoint,
      n         = as.integer(d$n),
      k         = as.integer(d[[ep$k_col]]),
      stringsAsFactors = FALSE
    )
  } else {
    d <- pick_primary_timepoint_continuous(raw, ep, ep$target_wk)
    if (!nrow(d)) return(data.frame())
    d <- attach_derived_change(d, raw, ep)
    d <- d[!is.na(d$change_mean) & !is.na(d$change_sd) &
           !is.na(d$n) & d$n > 1 & !is.na(d$drug) & nzchar(d$drug),
           , drop = FALSE]
    if (!nrow(d)) return(data.frame())
    data.frame(
      ref_id    = d$ref_id,
      trial     = d$trial,
      arm_no    = d$arm_no,
      drug      = d$drug,
      timepoint = d$timepoint,
      n         = as.integer(d$n),
      mean      = d$change_mean,
      sd        = d$change_sd,
      stringsAsFactors = FALSE
    )
  }
}

# ---------------------------------------------------------------------------
# Arm aggregation (a trial may carry several arms of the same drug, e.g. dose
# variants). Binary: sum n + k. Continuous: inverse-variance weighted mean and
# pooled SD, treating the variants as one arm.
# ---------------------------------------------------------------------------
aggregate_arms_binary <- function(rows) {
  if (nrow(rows) == 1L) return(list(n = rows$n[1], k = rows$k[1]))
  list(n = sum(rows$n), k = sum(rows$k))
}

aggregate_arms_continuous <- function(rows) {
  if (nrow(rows) == 1L)
    return(list(n = rows$n[1], mean = rows$mean[1], sd = rows$sd[1]))
  w  <- rows$n / (rows$sd^2)
  mu <- sum(w * rows$mean) / sum(w)
  ss <- sum((rows$n - 1) * rows$sd^2 + rows$n * (rows$mean - mu)^2)
  v  <- ss / (sum(rows$n) - 1)
  list(n = sum(rows$n), mean = mu, sd = sqrt(v))
}

aggregate_per_trial_drug <- function(arm_tbl, kind) {
  key <- paste(arm_tbl$ref_id, arm_tbl$drug, sep = "|")
  if (!anyDuplicated(key)) return(arm_tbl)
  out <- split(arm_tbl, key, drop = TRUE)
  collapsed <- lapply(out, function(rows) {
    if (nrow(rows) == 1L) return(rows)
    if (kind == "binary") {
      agg <- aggregate_arms_binary(rows)
      rows[1, "n"] <- agg$n; rows[1, "k"] <- agg$k
    } else {
      agg <- aggregate_arms_continuous(rows)
      rows[1, "n"] <- agg$n; rows[1, "mean"] <- agg$mean; rows[1, "sd"] <- agg$sd
    }
    rows[1, , drop = FALSE]
  })
  do.call(rbind, collapsed)
}

# ---------------------------------------------------------------------------
# prep_* : build the prepared data frame a model function consumes. These stop
# before any model is fit -- they are pure data shaping.
# ---------------------------------------------------------------------------

# Per-trial two-arm contrast input for one (drug_a, drug_b). NULL if no trial
# carries both drugs for this outcome.
prep_pairwise <- function(arm_tbl, drug_a, drug_b, kind) {
  if (!nrow(arm_tbl)) return(NULL)
  ta <- arm_tbl[arm_tbl$drug == drug_a, , drop = FALSE]
  tb <- arm_tbl[arm_tbl$drug == drug_b, , drop = FALSE]
  trials <- intersect(unique(ta$ref_id), unique(tb$ref_id))
  if (!length(trials)) return(NULL)
  per_trial <- vector("list", length(trials))
  for (i in seq_along(trials)) {
    rid <- trials[i]
    ra  <- ta[ta$ref_id == rid, , drop = FALSE]
    rb  <- tb[tb$ref_id == rid, , drop = FALSE]
    trial_name <- ra$trial[1] %||% rb$trial[1]
    if (kind == "binary") {
      a <- aggregate_arms_binary(ra); b <- aggregate_arms_binary(rb)
      per_trial[[i]] <- data.frame(
        ref_id = rid, trial = trial_name,
        event_a = a$k, n_a = a$n, event_b = b$k, n_b = b$n,
        stringsAsFactors = FALSE)
    } else {
      a <- aggregate_arms_continuous(ra); b <- aggregate_arms_continuous(rb)
      per_trial[[i]] <- data.frame(
        ref_id = rid, trial = trial_name,
        mean_a = a$mean, sd_a = a$sd, n_a = a$n,
        mean_b = b$mean, sd_b = b$sd, n_b = b$n,
        stringsAsFactors = FALSE)
    }
  }
  pt <- do.call(rbind, per_trial)
  if (!nrow(pt)) NULL else pt
}

# Per-trial k/n for a single drug's single-arm proportion. NULL if none.
prep_proportion <- function(arm_tbl, drug) {
  if (!nrow(arm_tbl)) return(NULL)
  d <- arm_tbl[arm_tbl$drug == drug & !is.na(arm_tbl$k), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  agg <- aggregate(d[, c("k", "n")],
                   by = list(ref_id = d$ref_id, trial = d$trial), FUN = sum)
  if (!nrow(agg)) NULL else agg
}

# Drug-aggregated arm-level long table for the network (>= 2 drugs per trial).
# NULL if no trial connects two treatments.
prep_nma <- function(arm_tbl, kind) {
  if (!nrow(arm_tbl)) return(NULL)
  arm_tbl <- aggregate_per_trial_drug(arm_tbl, kind)
  by_trial <- split(arm_tbl, arm_tbl$ref_id)
  ok_ids <- names(by_trial)[vapply(by_trial,
                                   function(d) length(unique(d$drug)) >= 2L,
                                   logical(1))]
  if (!length(ok_ids)) return(NULL)
  arm_tbl <- arm_tbl[as.character(arm_tbl$ref_id) %in% ok_ids, , drop = FALSE]
  if (!nrow(arm_tbl)) NULL else arm_tbl
}

# ---------------------------------------------------------------------------
# Drive the build.
# ---------------------------------------------------------------------------
all_pairwise        <- list()
all_pairwise_trials <- list()
all_proportion      <- list()
all_proportion_trials <- list()
all_nma             <- list()
all_nma_estimates   <- list()

# Placebo is forced into the drug_b slot so RR reads "drug vs placebo".
orient_pair <- function(p) if (identical(p[1], "Placebo")) c(p[2], p[1]) else p

cat(sprintf("Building meta-analysis tables (families: %s)...\n",
            paste(RUN_FAMILIES, collapse = ", ")))

for (spec in analyses) {
  # Build each outcome's arm table once; reused across families.
  arm_by_outcome <- setNames(
    lapply(spec$outcomes, function(oc) build_arm_table(spec_ep(spec, oc))),
    vapply(spec$outcomes, function(oc) oc$code, character(1)))
  codes <- names(arm_by_outcome)
  nonempty <- codes[vapply(arm_by_outcome, nrow, integer(1)) > 0]
  if (!length(nonempty)) {
    cat(sprintf("  %s: no contributing arms - skipped\n", spec$group))
    next
  }

  fn_pairwise   <- spec$model_pairwise   %||% model_pairwise
  fn_proportion <- spec$model_proportion %||% model_proportion
  fn_nma        <- spec$model_nma        %||% model_nma

  msg <- character(0)

  # --- pairwise -----------------------------------------------------------
  if ("pairwise" %in% RUN_FAMILIES) {
    pair_set <- list()
    for (code in nonempty) {
      by_trial <- split(arm_by_outcome[[code]]$drug, arm_by_outcome[[code]]$ref_id)
      for (ds in by_trial) {
        ds <- sort(unique(ds))
        if (length(ds) < 2) next
        cmb <- utils::combn(ds, 2)
        for (j in seq_len(ncol(cmb)))
          pair_set[[paste(cmb[1, j], cmb[2, j], sep = "|")]] <- c(cmb[1, j], cmb[2, j])
      }
    }
    pair_set <- lapply(pair_set, orient_pair)
    n_pw <- 0L
    for (pair in pair_set) {
      by_outcome <- setNames(
        lapply(codes, function(code)
          prep_pairwise(arm_by_outcome[[code]], pair[1], pair[2], spec$kind)),
        codes)
      if (all(vapply(by_outcome, is.null, logical(1)))) next
      res <- fn_pairwise(list(drug_a = pair[1], drug_b = pair[2],
                              by_outcome = by_outcome), spec)
      if (is.null(res) || !nrow(res$summary)) next
      s <- res$summary
      s$comparison_id <- paste(s$endpoint_group, s$outcome_code,
                               s$drug_a, s$drug_b, sep = "|")
      d <- res$detail
      cid <- paste(spec$group, d$outcome_code, pair[1], pair[2], sep = "|")
      d$outcome_code <- NULL
      d <- cbind(comparison_id = cid, d, stringsAsFactors = FALSE)
      all_pairwise[[length(all_pairwise) + 1L]] <- s
      all_pairwise_trials[[length(all_pairwise_trials) + 1L]] <- d
      n_pw <- n_pw + nrow(s)
    }
    msg <- c(msg, sprintf("%d pairwise", n_pw))
  }

  # --- single-arm proportion (binary only) --------------------------------
  if ("proportion" %in% RUN_FAMILIES && spec$kind == "binary") {
    drugs <- sort(unique(unlist(lapply(arm_by_outcome[nonempty], `[[`, "drug"))))
    n_prop <- 0L
    for (drug in drugs) {
      by_outcome <- setNames(
        lapply(codes, function(code) prep_proportion(arm_by_outcome[[code]], drug)),
        codes)
      if (all(vapply(by_outcome, is.null, logical(1)))) next
      res <- fn_proportion(list(drug = drug, by_outcome = by_outcome), spec)
      if (is.null(res) || !nrow(res$summary)) next
      s <- res$summary
      s$proportion_id <- paste(s$endpoint_group, s$outcome_code, s$drug, sep = "|")
      d <- res$detail
      pid <- paste(spec$group, d$outcome_code, drug, sep = "|")
      d$outcome_code <- NULL
      d <- cbind(proportion_id = pid, d, stringsAsFactors = FALSE)
      all_proportion[[length(all_proportion) + 1L]] <- s
      all_proportion_trials[[length(all_proportion_trials) + 1L]] <- d
      n_prop <- n_prop + nrow(s)
    }
    msg <- c(msg, sprintf("%d proportion", n_prop))
  }

  # --- network meta-analysis ----------------------------------------------
  if ("nma" %in% RUN_FAMILIES) {
    by_outcome <- setNames(
      lapply(codes, function(code) prep_nma(arm_by_outcome[[code]], spec$kind)),
      codes)
    res <- fn_nma(list(by_outcome = by_outcome), spec)
    n_net <- 0L; n_est <- 0L
    if (!is.null(res) && nrow(res$summary)) {
      s <- res$summary
      s$network_id <- paste(s$endpoint_group, s$outcome_code, sep = "|")
      all_nma[[length(all_nma) + 1L]] <- s
      n_net <- nrow(s)
      d <- res$detail
      if (!is.null(d) && nrow(d)) {
        nid <- paste(spec$group, d$outcome_code, sep = "|")
        d$outcome_code <- NULL
        d$network_id <- nid
        all_nma_estimates[[length(all_nma_estimates) + 1L]] <- d
        n_est <- nrow(d)
      }
    }
    msg <- c(msg, sprintf("%d networks, %d estimates", n_net, n_est))
  }

  cat(sprintf("  %s: %s\n", spec$group, paste(msg, collapse = "; ")))
}

# ---------------------------------------------------------------------------
# Assemble + write. Only families that ran are (over)written; the rest keep
# their existing rows (and built_at) untouched.
# ---------------------------------------------------------------------------
bind_or_empty <- function(lst, cols) {
  if (length(lst)) do.call(rbind, lst)
  else as.data.frame(setNames(replicate(length(cols), logical(0),
                                        simplify = FALSE), cols),
                     stringsAsFactors = FALSE)
}

built_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")
add_built_at <- function(df) { df$built_at <- if (nrow(df)) built_at else character(0); df }

# Final column orders (contract cols + driver id), mirroring the original tables.
COLS_PAIRWISE        <- c(.PAIRWISE_SUMMARY_COLS, "comparison_id")
COLS_PAIRWISE_TRIALS <- c("comparison_id", setdiff(.PAIRWISE_DETAIL_COLS, "outcome_code"))
COLS_PROP            <- c(.PROP_SUMMARY_COLS, "proportion_id")
COLS_PROP_TRIALS     <- c("proportion_id", setdiff(.PROP_DETAIL_COLS, "outcome_code"))
COLS_NMA             <- c(.NMA_SUMMARY_COLS, "network_id")
COLS_NMA_ESTIMATES   <- c(setdiff(.NMA_DETAIL_COLS, "outcome_code"), "network_id")

# table name -> (family, accumulated rows, column template, index spec)
writes <- list(
  ma_pairwise = list(family = "pairwise", rows = all_pairwise, cols = COLS_PAIRWISE,
    idx = "CREATE INDEX IF NOT EXISTS idx_ma_pairwise_lookup
           ON ma_pairwise(endpoint_group, outcome_code, drug_a, drug_b)"),
  ma_pairwise_trials = list(family = "pairwise", rows = all_pairwise_trials,
    cols = COLS_PAIRWISE_TRIALS,
    idx = "CREATE INDEX IF NOT EXISTS idx_ma_pairwise_trials_cid
           ON ma_pairwise_trials(comparison_id)"),
  ma_proportion = list(family = "proportion", rows = all_proportion, cols = COLS_PROP,
    idx = "CREATE INDEX IF NOT EXISTS idx_ma_proportion_lookup
           ON ma_proportion(endpoint_group, outcome_code, drug)"),
  ma_proportion_trials = list(family = "proportion", rows = all_proportion_trials,
    cols = COLS_PROP_TRIALS,
    idx = "CREATE INDEX IF NOT EXISTS idx_ma_proportion_trials_pid
           ON ma_proportion_trials(proportion_id)"),
  ma_nma = list(family = "nma", rows = all_nma, cols = COLS_NMA,
    idx = "CREATE INDEX IF NOT EXISTS idx_ma_nma_lookup
           ON ma_nma(endpoint_group, outcome_code)"),
  ma_nma_estimates = list(family = "nma", rows = all_nma_estimates,
    cols = COLS_NMA_ESTIMATES,
    idx = "CREATE INDEX IF NOT EXISTS idx_ma_nma_estimates_lookup
           ON ma_nma_estimates(network_id, drug_a, drug_b)")
)

cat(sprintf("\nWriting tables to %s ...\n", DB_PATH))
written <- character(0)
for (tbl in names(writes)) {
  w <- writes[[tbl]]
  if (!(w$family %in% RUN_FAMILIES)) {
    cat(sprintf("  %-22s skipped (family '%s' not selected)\n", tbl, w$family))
    next
  }
  df <- add_built_at(bind_or_empty(w$rows, w$cols))
  dbWriteTable(con, tbl, df, overwrite = TRUE)
  dbExecute(con, w$idx)
  written <- c(written, sprintf("%s: %d", tbl, nrow(df)))
}

cat(sprintf("Done. %s\n", paste(written, collapse = ", ")))
