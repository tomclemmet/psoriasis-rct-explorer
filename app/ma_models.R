# ===========================================================================
# ma_models.R  --  THE FILE YOU EDIT TO CHANGE / ADD STATISTICAL MODELS.
#
# meta_analyse.R (the driver) does all the plumbing: it reads the views, picks
# primary timepoints, builds per-arm tables, enumerates comparisons, and writes
# the six ma_* tables the Shiny app reads. It hands each "model" function a
# clean, prepared data frame and asks for results back in a fixed structure.
# You own the model functions below; swap in JAGS/STAN/brms whenever you like.
#
# ---------------------------------------------------------------------------
# WHAT A MODEL FUNCTION RECEIVES
# ---------------------------------------------------------------------------
# Each family's model fn is called once per analysis *unit* with:
#   pairwise:    unit = list(drug_a, drug_b, by_outcome = <named list>)
#   proportion:  unit = list(drug,           by_outcome = <named list>)
#   nma:         unit = list(               by_outcome = <named list>)
# and the `spec` it belongs to (group / view / kind / target_wk / outcomes).
#
# `by_outcome` is keyed by outcome_code; each element is the prepared data frame
# for that one outcome (built by the driver's prep_* functions):
#   pairwise binary:     per-trial event_a/n_a/event_b/n_b
#   pairwise continuous: per-trial mean_a/sd_a/n_a, mean_b/sd_b/n_b
#   proportion:          per-trial k/n (one drug)
#   nma:                 drug-aggregated arm-level long table (drug,n,k|mean,sd)
#
# The default models loop `by_outcome` and fit each outcome independently -- the
# *exact* frequentist analysis the app shipped with. A joint model (e.g. an
# ordered-multinomial PASI model) ignores the loop, fits all outcomes at once,
# and returns rows for every threshold. Nothing downstream changes: the driver
# keys the SQLite ids off the `outcome_code` your rows carry.
#
# ---------------------------------------------------------------------------
# WHAT A MODEL FUNCTION RETURNS  (the contract)
# ---------------------------------------------------------------------------
# A result built with pairwise_result() / proportion_result() / nma_result():
#   list(summary = <data frame, 1 row per outcome>,
#        detail  = <data frame, the forest / per-study rows>)
# Every row carries `outcome_code` so one result may span several outcomes.
# Missing OPTIONAL columns are filled with NA; missing ESSENTIAL columns are a
# hard error (so a malformed hand-built Bayesian result fails loudly, early).
# Return NULL to contribute nothing (no contributing data for this unit).
#
# ---------------------------------------------------------------------------
# THREE SHARP EDGES THE APP SILENTLY ASSUMES -- honour them or plots go wrong
# ---------------------------------------------------------------------------
# 1. SCALE. Pairwise/NMA relative-risk estimates are stored on the LOG scale
#    (app.R back-transforms with exp()). Mean differences are stored as-is.
#    Single-arm proportions are stored already back-transformed to the natural
#    0-1 scale (apply inv_logit() before returning). model_scale(spec) tells you
#    which: "rr" | "md" | "prop".
# 2. FE vs RE. The app toggles between the *_fe and *_re columns. A Bayesian fit
#    with a single posterior estimate should put it in BOTH -- use both(te,lo,hi)
#    which returns all six columns -- or one toggle state shows blanks.
#
# (Forest-plot square sizes are NOT a model concern: the app sizes per-trial
# squares from the trial's sample size n, and draws network/pooled rows uniform.
# Detail rows just need to carry n -- no study weights are stored or required.)
# ===========================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

inv_logit <- function(x) exp(x) / (1 + exp(x))

# Which storage scale this spec's estimates use (see sharp edge #1).
model_scale <- function(spec) {
  if (identical(spec$family, "proportion")) return("prop")
  if (identical(spec$kind, "binary")) "rr" else "md"
}

# Put a single estimate into both the fixed- and random-effect slots.
both <- function(te, lo, hi) {
  list(te_fe = te, lo_fe = lo, hi_fe = hi,
       te_re = te, lo_re = lo, hi_re = hi)
}

# ---------------------------------------------------------------------------
# Result contract: canonical column orders (mirror the existing ma_* tables so
# app.R needs no changes) + constructors that validate essentials and fill the
# rest with NA. The driver appends the id (comparison_id / proportion_id /
# network_id) and built_at; detail frames carry outcome_code as a join key that
# the driver consumes then drops.
# ---------------------------------------------------------------------------
.PAIRWISE_SUMMARY_COLS <- c(
  "drug_a", "drug_b", "endpoint_group", "outcome_code", "sm", "n_studies",
  "te_fe", "se_fe", "lo_fe", "hi_fe", "z_fe", "p_fe",
  "te_re", "se_re", "lo_re", "hi_re", "z_re", "p_re",
  "tau2", "i2", "q", "q_df", "q_pval", "method_tau")
.PAIRWISE_DETAIL_COLS <- c(
  "ref_id", "trial", "te", "se", "lo", "hi",
  "event_a", "n_a", "event_b", "n_b", "mean_a", "sd_a", "mean_b", "sd_b",
  "outcome_code")

.PROP_SUMMARY_COLS <- c(
  "drug", "endpoint_group", "outcome_code", "sm", "n_studies",
  "te_fe", "lo_fe", "hi_fe", "te_re", "lo_re", "hi_re",
  "tau2", "i2", "q", "q_df", "q_pval", "method_tau")
.PROP_DETAIL_COLS <- c(
  "ref_id", "trial", "k", "n", "p", "lo", "hi",
  "outcome_code")

.NMA_SUMMARY_COLS <- c(
  "endpoint_group", "outcome_code", "sm", "status", "n_studies",
  "n_treatments", "n_pairwise", "tau2", "i2",
  "q_total", "q_het", "q_inc", "p_inc")
.NMA_DETAIL_COLS <- c(
  "drug_a", "drug_b", "te_fe", "se_fe", "lo_fe", "hi_fe",
  "te_re", "se_re", "lo_re", "hi_re", "n_direct", "n_indirect",
  "outcome_code")

# Coerce a data frame to a contract: error on missing essentials, NA-fill the
# rest, return columns in canonical order.
.fit_contract <- function(df, cols, essential, label) {
  if (is.null(df) || !nrow(df)) {
    empty <- as.data.frame(setNames(
      replicate(length(cols), logical(0), simplify = FALSE), cols),
      stringsAsFactors = FALSE)
    return(empty)
  }
  missing_ess <- setdiff(essential, names(df))
  if (length(missing_ess))
    stop(sprintf("%s missing required column(s): %s",
                 label, paste(missing_ess, collapse = ", ")), call. = FALSE)
  for (col in setdiff(cols, names(df))) df[[col]] <- NA
  df[, cols, drop = FALSE]
}

pairwise_result <- function(summary, detail) {
  list(
    summary = .fit_contract(summary, .PAIRWISE_SUMMARY_COLS,
      c("drug_a", "drug_b", "endpoint_group", "outcome_code",
        "te_fe", "lo_fe", "hi_fe", "te_re", "lo_re", "hi_re"),
      "pairwise summary"),
    detail = .fit_contract(detail, .PAIRWISE_DETAIL_COLS,
      c("outcome_code", "trial", "te", "lo", "hi"),
      "pairwise detail"))
}

proportion_result <- function(summary, detail) {
  list(
    summary = .fit_contract(summary, .PROP_SUMMARY_COLS,
      c("drug", "endpoint_group", "outcome_code",
        "te_fe", "lo_fe", "te_re", "lo_re", "hi_re"),
      "proportion summary"),
    detail = .fit_contract(detail, .PROP_DETAIL_COLS,
      c("outcome_code", "trial", "k", "n", "p", "lo", "hi"),
      "proportion detail"))
}

nma_result <- function(summary, detail) {
  list(
    summary = .fit_contract(summary, .NMA_SUMMARY_COLS,
      c("endpoint_group", "outcome_code", "status", "i2", "n_studies"),
      "nma summary"),
    detail = .fit_contract(detail, .NMA_DETAIL_COLS,
      c("outcome_code", "drug_a", "drug_b",
        "te_fe", "lo_fe", "hi_fe", "te_re", "lo_re", "hi_re", "n_direct"),
      "nma detail"))
}

# ---------------------------------------------------------------------------
# tidy_* helpers: turn a fitted meta::/netmeta:: object into contract rows.
# These are the one-liners the default models below lean on; copy them as a
# template when wiring up a non-meta model (build the same columns by hand).
# ---------------------------------------------------------------------------
tidy_pairwise <- function(ma, pt, drug_a, drug_b, group, outcome, kind) {
  summary <- data.frame(
    drug_a = drug_a, drug_b = drug_b,
    endpoint_group = group, outcome_code = outcome,
    sm = ma$sm, n_studies = ma$k,
    te_fe = ma$TE.common, se_fe = ma$seTE.common,
    lo_fe = ma$lower.common, hi_fe = ma$upper.common,
    z_fe = ma$zval.common, p_fe = ma$pval.common,
    te_re = ma$TE.random, se_re = ma$seTE.random,
    lo_re = ma$lower.random, hi_re = ma$upper.random,
    z_re = ma$zval.random, p_re = ma$pval.random,
    tau2 = ma$tau2 %||% NA_real_, i2 = ma$I2,
    q = ma$Q, q_df = ma$df.Q, q_pval = ma$pval.Q,
    method_tau = "REML", stringsAsFactors = FALSE)

  detail <- data.frame(
    ref_id = pt$ref_id, trial = pt$trial,
    te = ma$TE, se = ma$seTE, lo = ma$lower, hi = ma$upper,
    stringsAsFactors = FALSE)
  if (kind == "binary") {
    detail$event_a <- pt$event_a; detail$n_a <- pt$n_a
    detail$event_b <- pt$event_b; detail$n_b <- pt$n_b
    detail$mean_a <- NA_real_; detail$sd_a <- NA_real_
    detail$mean_b <- NA_real_; detail$sd_b <- NA_real_
  } else {
    detail$event_a <- NA_integer_; detail$n_a <- pt$n_a
    detail$event_b <- NA_integer_; detail$n_b <- pt$n_b
    detail$mean_a <- pt$mean_a; detail$sd_a <- pt$sd_a
    detail$mean_b <- pt$mean_b; detail$sd_b <- pt$sd_b
  }
  detail$outcome_code <- outcome
  list(summary = summary, detail = detail)
}

tidy_proportion <- function(ma, agg, drug, group, outcome) {
  summary <- data.frame(
    drug = drug, endpoint_group = group, outcome_code = outcome,
    sm = "proportion", n_studies = ma$k,
    te_fe = inv_logit(ma$TE.common),
    lo_fe = inv_logit(ma$lower.common), hi_fe = inv_logit(ma$upper.common),
    te_re = inv_logit(ma$TE.random),
    lo_re = inv_logit(ma$lower.random), hi_re = inv_logit(ma$upper.random),
    tau2 = ma$tau2 %||% NA_real_, i2 = ma$I2,
    q = ma$Q, q_df = ma$df.Q, q_pval = ma$pval.Q,
    method_tau = "REML", stringsAsFactors = FALSE)

  p_hat <- agg$k / agg$n
  z <- 1.959964
  denom  <- 1 + z^2 / agg$n
  centre <- (p_hat + z^2 / (2 * agg$n)) / denom
  half   <- z * sqrt(p_hat * (1 - p_hat) / agg$n + z^2 / (4 * agg$n^2)) / denom
  detail <- data.frame(
    ref_id = agg$ref_id, trial = agg$trial, k = agg$k, n = agg$n, p = p_hat,
    lo = pmax(0, centre - half), hi = pmin(1, centre + half),
    outcome_code = outcome, stringsAsFactors = FALSE)
  list(summary = summary, detail = detail)
}

tidy_netmeta <- function(nm, group, outcome) {
  summary <- data.frame(
    endpoint_group = group, outcome_code = outcome,
    sm = nm$sm, status = "ok",
    n_studies = nm$k, n_treatments = length(nm$trts), n_pairwise = nm$m,
    tau2 = nm$tau2 %||% NA_real_, i2 = nm$I2 %||% NA_real_,
    q_total = nm$Q %||% NA_real_, q_het = nm$Q.heterogeneity %||% NA_real_,
    q_inc = nm$Q.inconsistency %||% NA_real_,
    p_inc = nm$pval.Q.inconsistency %||% NA_real_,
    stringsAsFactors = FALSE)

  treats <- nm$trts
  n_t <- length(treats)
  pairs <- expand.grid(a = seq_len(n_t), b = seq_len(n_t),
                       KEEP.OUT.ATTRS = FALSE)
  pairs <- pairs[pairs$a != pairs$b, , drop = FALSE]
  pw <- nm$data
  count_direct <- function(ta, tb)
    sum((pw$.treat1 == ta & pw$.treat2 == tb) |
        (pw$.treat1 == tb & pw$.treat2 == ta))
  detail <- data.frame(
    drug_a = treats[pairs$a], drug_b = treats[pairs$b],
    te_fe = nm$TE.common[cbind(pairs$a, pairs$b)],
    se_fe = nm$seTE.common[cbind(pairs$a, pairs$b)],
    lo_fe = nm$lower.common[cbind(pairs$a, pairs$b)],
    hi_fe = nm$upper.common[cbind(pairs$a, pairs$b)],
    te_re = nm$TE.random[cbind(pairs$a, pairs$b)],
    se_re = nm$seTE.random[cbind(pairs$a, pairs$b)],
    lo_re = nm$lower.random[cbind(pairs$a, pairs$b)],
    hi_re = nm$upper.random[cbind(pairs$a, pairs$b)],
    stringsAsFactors = FALSE)
  detail$n_direct   <- mapply(count_direct, detail$drug_a, detail$drug_b)
  detail$n_indirect <- pmax(0L, nm$k - detail$n_direct)
  detail$outcome_code <- outcome
  list(summary = summary, detail = detail)
}

# Skeleton summary row for a network too sparse / disconnected to fit.
tidy_nma_sparse <- function(group, outcome, sm) {
  data.frame(
    endpoint_group = group, outcome_code = outcome,
    sm = sm, status = "sparse",
    n_studies = NA_integer_, n_treatments = NA_integer_,
    n_pairwise = NA_integer_, tau2 = NA_real_, i2 = NA_real_,
    q_total = NA_real_, q_het = NA_real_, q_inc = NA_real_, p_inc = NA_real_,
    stringsAsFactors = FALSE)
}

# Helpers used only by the default NMA model: contrast-level long format from
# an arm-level table, and the netmeta fit + connectivity guard.
.nma_pairwise <- function(arm, kind) {
  args <- list(treat = arm$drug, studlab = arm$ref_id,
               data = arm, allstudies = TRUE)
  if (kind == "binary")
    args <- c(args, list(event = arm$k, n = arm$n, sm = "RR"))
  else
    args <- c(args, list(n = arm$n, mean = arm$mean, sd = arm$sd, sm = "MD"))
  pw <- tryCatch(suppressWarnings(do.call(meta::pairwise, args)),
                 error = function(e) NULL)
  if (is.null(pw) || !nrow(pw)) return(NULL)
  pw <- pw[is.finite(pw$TE) & is.finite(pw$seTE) & pw$seTE > 0, , drop = FALSE]
  if (!nrow(pw)) return(NULL)
  pw
}

.nma_fit <- function(pw) {
  treats <- unique(c(pw$treat1, pw$treat2))
  if (length(treats) < 2L) return(NULL)
  nm <- tryCatch(suppressWarnings(netmeta::netmeta(
    TE = pw$TE, seTE = pw$seTE,
    treat1 = pw$treat1, treat2 = pw$treat2,
    studlab = pw$studlab, data = pw,
    sm = pw$sm[1] %||% attr(pw, "sm") %||% "RR",
    common = TRUE, random = TRUE,
    reference.group = if ("Placebo" %in% treats) "Placebo" else treats[1],
    tol.multiarm = 0.5)), error = function(e) NULL)
  if (is.null(nm)) return(NULL)
  if (!isTRUE(nm$n.subnets == 1L) && !is.null(nm$n.subnets)) return(NULL)
  nm
}

# ===========================================================================
# DEFAULT MODELS  --  reproduce the app's original frequentist analyses.
# Each loops the spec's outcomes and fits them independently. Replace the body
# of any of these (or set spec$model_<family> to your own fn) to take over.
# ===========================================================================
model_pairwise <- function(unit, spec) {
  summaries <- list(); details <- list()
  for (code in names(unit$by_outcome)) {
    pt <- unit$by_outcome[[code]]
    if (is.null(pt) || !nrow(pt)) next
    ma <- tryCatch(suppressWarnings(
      if (spec$kind == "binary")
        meta::metabin(event.e = pt$event_a, n.e = pt$n_a,
                      event.c = pt$event_b, n.c = pt$n_b,
                      studlab = pt$trial, sm = "RR", method = "Inverse",
                      method.tau = "REML", common = TRUE, random = TRUE,
                      warn = FALSE)
      else
        meta::metacont(n.e = pt$n_a, mean.e = pt$mean_a, sd.e = pt$sd_a,
                       n.c = pt$n_b, mean.c = pt$mean_b, sd.c = pt$sd_b,
                       studlab = pt$trial, sm = "MD", method.tau = "REML",
                       common = TRUE, random = TRUE, warn = FALSE)
    ), error = function(e) NULL)
    if (is.null(ma)) next
    tid <- tidy_pairwise(ma, pt, unit$drug_a, unit$drug_b,
                         spec$group, code, spec$kind)
    summaries[[code]] <- tid$summary
    details[[code]]   <- tid$detail
  }
  if (!length(summaries)) return(NULL)
  pairwise_result(do.call(rbind, summaries), do.call(rbind, details))
}

model_proportion <- function(unit, spec) {
  summaries <- list(); details <- list()
  for (code in names(unit$by_outcome)) {
    agg <- unit$by_outcome[[code]]
    if (is.null(agg) || !nrow(agg)) next
    ma <- tryCatch(suppressWarnings(
      meta::metaprop(event = agg$k, n = agg$n, studlab = agg$trial,
                     sm = "PLOGIT", method = "Inverse", method.tau = "REML",
                     common = TRUE, random = TRUE, warn = FALSE)),
      error = function(e) NULL)
    if (is.null(ma)) next
    tid <- tidy_proportion(ma, agg, unit$drug, spec$group, code)
    summaries[[code]] <- tid$summary
    details[[code]]   <- tid$detail
  }
  if (!length(summaries)) return(NULL)
  proportion_result(do.call(rbind, summaries), do.call(rbind, details))
}

model_nma <- function(unit, spec) {
  summaries <- list(); details <- list()
  sm <- if (spec$kind == "binary") "RR" else "MD"
  for (code in names(unit$by_outcome)) {
    arm <- unit$by_outcome[[code]]
    pw  <- if (is.null(arm) || !nrow(arm)) NULL else .nma_pairwise(arm, spec$kind)
    nm  <- if (is.null(pw)) NULL else .nma_fit(pw)
    if (is.null(nm)) {
      summaries[[code]] <- tidy_nma_sparse(spec$group, code, sm)
      next
    }
    tid <- tidy_netmeta(nm, spec$group, code)
    summaries[[code]] <- tid$summary
    details[[code]]   <- tid$detail
  }
  if (!length(summaries)) return(NULL)
  det <- if (length(details)) do.call(rbind, details) else NULL
  nma_result(do.call(rbind, summaries), det)
}
