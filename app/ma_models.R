# ===========================================================================
# ma_models.R  --  THE FILE YOU EDIT TO CHANGE / ADD MODELS.
#
# The driver (meta_analyse.R) hands each model a prepared data frame and stores
# whatever it returns. Each model does two things: (1) FIT, (2) POST-PROCESS the
# fit into the result format. Swap in JAGS/STAN/brms freely -- only the returned
# columns matter. The format + validators live in ma_contract.R (you rarely
# touch that); the *_result() constructors below come from there.
#
# ---------------------------------------------------------------------------
# WHAT YOU RECEIVE
# ---------------------------------------------------------------------------
#   model_pairwise(unit, spec)    unit = list(drug_a, drug_b, by_outcome)
#   model_proportion(unit, spec)  unit = list(drug,           by_outcome)
#   model_nma(unit, spec)         unit = list(               by_outcome)
# `spec` carries group / view / kind ("binary"|"continuous") / target_wk.
# `by_outcome` is a named list keyed by outcome_code; each element is the
# prepared per-outcome data frame:
#   pairwise binary      per-trial event_a/n_a, event_b/n_b
#   pairwise continuous  per-trial mean_a/sd_a/n_a, mean_b/sd_b/n_b
#   proportion           per-trial k/n (a single drug)
#   nma                  drug-aggregated arm-level rows (drug, n, k|mean,sd)
# Looping by_outcome and fitting each outcome on its own reproduces today's
# analysis; a joint model (e.g. ordinal PASI) can fit them together and emit
# rows for every threshold in one go -- the driver keys ids off outcome_code.
#
# ---------------------------------------------------------------------------
# WHAT YOU RETURN
# ---------------------------------------------------------------------------
#   pairwise_result(summary, detail) | proportion_result(..) | nma_result(..)
# summary = one row per outcome, detail = the per-study / forest rows; both
# data frames carry outcome_code. Columns are listed in ma_contract.R: missing
# optional ones NA-fill, a missing ESSENTIAL one is a hard error. Return NULL to
# contribute nothing for this unit.
#
# TWO THINGS THE APP ASSUMES (or the plots come out wrong):
#   * SCALE: store relative risks on the LOG scale (the app exp()s them); mean
#     differences as-is; single-arm proportions back-transformed to 0-1 with
#     inv_logit(). model_scale(spec) returns "rr" | "md" | "prop".
#   * FE & RE: the app toggles the *_fe vs *_re columns -- fill BOTH. A single
#     Bayesian posterior estimate goes in both slots; both(te, lo, hi) returns
#     the six columns for you.
# ===========================================================================

`%||%`    <- function(a, b) if (is.null(a)) b else a
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Put one estimate into both the fixed- and random-effect slots.
both <- function(te, lo, hi)
  list(te_fe = te, lo_fe = lo, hi_fe = hi, te_re = te, lo_re = lo, hi_re = hi)

# Storage scale for this spec (see SCALE above).
model_scale <- function(spec) {
  if (identical(spec$family, "proportion")) return("prop")
  if (identical(spec$kind, "binary")) "rr" else "md"
}

# ===========================================================================
# PAIRWISE  --  one meta-analysis per drug-pair per outcome.
# Binary -> metabin (RR, log scale); continuous -> metacont (MD).
# ===========================================================================
model_pairwise <- function(unit, spec) {
  summaries <- list(); details <- list()
  for (code in names(unit$by_outcome)) {
    pt <- unit$by_outcome[[code]]
    if (is.null(pt) || !nrow(pt)) next

    ## ---- FIT -------------------------------------------------------------
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

    ## ---- POST-PROCESS (RR stays on the log scale) ------------------------
    summaries[[code]] <- data.frame(
      drug_a = unit$drug_a, drug_b = unit$drug_b,
      endpoint_group = spec$group, outcome_code = code,
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

    det <- data.frame(
      ref_id = pt$ref_id, trial = pt$trial,
      te = ma$TE, se = ma$seTE, lo = ma$lower, hi = ma$upper,
      stringsAsFactors = FALSE)
    if (spec$kind == "binary") {
      det$event_a <- pt$event_a; det$n_a <- pt$n_a
      det$event_b <- pt$event_b; det$n_b <- pt$n_b
      det$mean_a <- NA_real_; det$sd_a <- NA_real_
      det$mean_b <- NA_real_; det$sd_b <- NA_real_
    } else {
      det$event_a <- NA_integer_; det$n_a <- pt$n_a
      det$event_b <- NA_integer_; det$n_b <- pt$n_b
      det$mean_a <- pt$mean_a; det$sd_a <- pt$sd_a
      det$mean_b <- pt$mean_b; det$sd_b <- pt$sd_b
    }
    det$outcome_code <- code
    details[[code]] <- det
  }
  if (!length(summaries)) return(NULL)
  pairwise_result(do.call(rbind, summaries), do.call(rbind, details))
}

# ===========================================================================
# PROPORTION  --  single-arm response rate per drug per outcome (metaprop).
# Pooled estimate stored back-transformed to 0-1; per-study Wilson CIs.
# ===========================================================================
model_proportion <- function(unit, spec) {
  summaries <- list(); details <- list()
  for (code in names(unit$by_outcome)) {
    agg <- unit$by_outcome[[code]]
    if (is.null(agg) || !nrow(agg)) next

    ## ---- FIT -------------------------------------------------------------
    ma <- tryCatch(suppressWarnings(
      meta::metaprop(event = agg$k, n = agg$n, studlab = agg$trial,
                     sm = "PLOGIT", method = "Inverse", method.tau = "REML",
                     common = TRUE, random = TRUE, warn = FALSE)),
      error = function(e) NULL)
    if (is.null(ma)) next

    ## ---- POST-PROCESS (back-transform pooled logit -> 0-1) ---------------
    summaries[[code]] <- data.frame(
      drug = unit$drug, endpoint_group = spec$group, outcome_code = code,
      sm = "proportion", n_studies = ma$k,
      te_fe = inv_logit(ma$TE.common),
      lo_fe = inv_logit(ma$lower.common), hi_fe = inv_logit(ma$upper.common),
      te_re = inv_logit(ma$TE.random),
      lo_re = inv_logit(ma$lower.random), hi_re = inv_logit(ma$upper.random),
      tau2 = ma$tau2 %||% NA_real_, i2 = ma$I2,
      q = ma$Q, q_df = ma$df.Q, q_pval = ma$pval.Q,
      method_tau = "REML", stringsAsFactors = FALSE)

    # Per-study point + Wilson score interval (independent of the pooling).
    p_hat  <- agg$k / agg$n
    z      <- 1.959964
    denom  <- 1 + z^2 / agg$n
    centre <- (p_hat + z^2 / (2 * agg$n)) / denom
    half   <- z * sqrt(p_hat * (1 - p_hat) / agg$n + z^2 / (4 * agg$n^2)) / denom
    details[[code]] <- data.frame(
      ref_id = agg$ref_id, trial = agg$trial, k = agg$k, n = agg$n, p = p_hat,
      lo = pmax(0, centre - half), hi = pmin(1, centre + half),
      outcome_code = code, stringsAsFactors = FALSE)
  }
  if (!length(summaries)) return(NULL)
  proportion_result(do.call(rbind, summaries), do.call(rbind, details))
}

# ===========================================================================
# NMA  --  one frequentist network meta-analysis per outcome (netmeta).
# The two helpers below prepare + fit the network; a custom (e.g. arm-based
# Bayesian) model would replace this whole section and fit the arm-level table
# unit$by_outcome[[code]] directly.
# ===========================================================================

# Arm-level table -> contrast-level long format (meta::pairwise), dropping
# non-finite contrasts. Returns NULL if nothing usable.
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

# Fit netmeta; return NULL unless the network is a single connected component
# (the app shows a "too sparse" note for those).
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

model_nma <- function(unit, spec) {
  summaries <- list(); details <- list()
  sm <- if (spec$kind == "binary") "RR" else "MD"
  for (code in names(unit$by_outcome)) {
    arm <- unit$by_outcome[[code]]

    ## ---- FIT -------------------------------------------------------------
    pw <- if (is.null(arm) || !nrow(arm)) NULL else .nma_pairwise(arm, spec$kind)
    nm <- if (is.null(pw)) NULL else .nma_fit(pw)

    ## ---- POST-PROCESS ----------------------------------------------------
    if (is.null(nm)) {                 # network too sparse / disconnected
      summaries[[code]] <- data.frame(
        endpoint_group = spec$group, outcome_code = code,
        sm = sm, status = "sparse",
        n_studies = NA_integer_, n_treatments = NA_integer_,
        n_pairwise = NA_integer_, tau2 = NA_real_, i2 = NA_real_,
        q_total = NA_real_, q_het = NA_real_, q_inc = NA_real_, p_inc = NA_real_,
        stringsAsFactors = FALSE)
      next
    }
    summaries[[code]] <- data.frame(
      endpoint_group = spec$group, outcome_code = code,
      sm = nm$sm, status = "ok",
      n_studies = nm$k, n_treatments = length(nm$trts), n_pairwise = nm$m,
      tau2 = nm$tau2 %||% NA_real_, i2 = nm$I2 %||% NA_real_,
      q_total = nm$Q %||% NA_real_, q_het = nm$Q.heterogeneity %||% NA_real_,
      q_inc = nm$Q.inconsistency %||% NA_real_,
      p_inc = nm$pval.Q.inconsistency %||% NA_real_,
      stringsAsFactors = FALSE)

    # One row per ordered drug pair, read off netmeta's effect matrices.
    treats <- nm$trts; n_t <- length(treats)
    pairs  <- expand.grid(a = seq_len(n_t), b = seq_len(n_t),
                          KEEP.OUT.ATTRS = FALSE)
    pairs  <- pairs[pairs$a != pairs$b, , drop = FALSE]
    pwd    <- nm$data
    count_direct <- function(ta, tb)
      sum((pwd$.treat1 == ta & pwd$.treat2 == tb) |
          (pwd$.treat1 == tb & pwd$.treat2 == ta))
    det <- data.frame(
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
    det$n_direct   <- mapply(count_direct, det$drug_a, det$drug_b)
    det$n_indirect <- pmax(0L, nm$k - det$n_direct)
    det$outcome_code <- code
    details[[code]] <- det
  }
  if (!length(summaries)) return(NULL)
  det <- if (length(details)) do.call(rbind, details) else NULL
  nma_result(do.call(rbind, summaries), det)
}
