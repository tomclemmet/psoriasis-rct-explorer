# Offline meta-analysis builder. Reads psoriasis-rcts.sqlite, computes
# every pooled estimate the Shiny app might display, writes the results
# back into the same SQLite file in four new tables:
#
#   ma_pairwise          - one row per (drug_a, drug_b, endpoint, outcome)
#                          comparison. FE + RE pooled estimate, plus
#                          heterogeneity stats.
#   ma_pairwise_trials   - per-trial contributions to each ma_pairwise row
#                          (the squares + tooltip payload of a forest plot).
#   ma_proportion        - one row per (drug, endpoint, outcome) pooled
#                          response rate (single-arm meta-analysis).
#   ma_proportion_trials - per-arm contributions to each ma_proportion row.
#
# Pairwise comparisons are computed for every (drug_a, drug_b) head-to-head
# that has at least one trial in the dataset, AND every (drug, Placebo)
# pair (so the no-filter view in the app can draw a forest plot of every
# drug vs placebo without recomputing).
#
# Per-trial primary timepoint: the latest measurement at or before a target
# week per (trial, arm, outcome). Targets: 16 wk for PASI, 24 wk for DLQI
# and safety. Trials reporting only later timepoints are excluded; trials
# reporting earlier-only timepoints contribute their latest available row.
#
# Effect measures:
#   binary endpoints (PASI thresholds, DLQI binaries, safety): RR
#   continuous endpoints (abs PASI, abs DLQI): MD on change-from-baseline
#   single-drug proportions: logit-transformed proportion (PLOGIT)
#
# Both common-effect (FE, inverse-variance) and random-effects (REML)
# estimates are stored on every row. The app picks which to emphasise.
#
# Re-running this script overwrites all four tables. Cheap (~hundreds of
# meta::* calls).

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

# meta:::settings.meta() chatter is just noise here.
suppressMessages(settings.meta(CIbracket = "(", print.tau2 = FALSE))

# ---------------------------------------------------------------------------
# Endpoint catalogue. Mirrors app/app.R's endpoint_groups but expressed in
# the column names of the v_* views directly. Each entry says: which view to
# read, which outcome columns to pool, the effect type (binary or
# continuous), and the per-trial primary-timepoint target in weeks.
# ---------------------------------------------------------------------------
endpoints <- list(
  list(group = "pasi", outcome = "pasi50",  view = "v_pasi", kind = "binary",
       k_col = "pasi50",  target_wk = 16),
  list(group = "pasi", outcome = "pasi75",  view = "v_pasi", kind = "binary",
       k_col = "pasi75",  target_wk = 16),
  list(group = "pasi", outcome = "pasi90",  view = "v_pasi", kind = "binary",
       k_col = "pasi90",  target_wk = 16),
  list(group = "pasi", outcome = "pasi100", view = "v_pasi", kind = "binary",
       k_col = "pasi100", target_wk = 16),
  list(group = "pasi_abs", outcome = "abs_pasi_change", view = "v_pasi_abs",
       kind = "continuous", mean_col = "abs_pasi_change_mean",
       sd_col = "abs_pasi_change_sd",
       baseline_mean_col = "abs_pasi_mean", baseline_sd_col = "abs_pasi_sd",
       followup_mean_col = "abs_pasi_mean", followup_sd_col = "abs_pasi_sd",
       target_wk = 16),
  list(group = "dlqi", outcome = "dlqi_0_1",     view = "v_dlqi", kind = "binary",
       k_col = "dlqi_0_1",     target_wk = 24),
  list(group = "dlqi", outcome = "dlqi_0",       view = "v_dlqi", kind = "binary",
       k_col = "dlqi_0",       target_wk = 24),
  list(group = "dlqi", outcome = "dlqi_le5",     view = "v_dlqi", kind = "binary",
       k_col = "dlqi_le5",     target_wk = 24),
  list(group = "dlqi", outcome = "dlqi_5pt_dec", view = "v_dlqi", kind = "binary",
       k_col = "dlqi_5pt_dec", target_wk = 24),
  list(group = "dlqi", outcome = "dlqi_4pt_dec", view = "v_dlqi", kind = "binary",
       k_col = "dlqi_4pt_dec", target_wk = 24),
  list(group = "dlqi_abs", outcome = "abs_dlqi_change", view = "v_dlqi",
       kind = "continuous", mean_col = "abs_dlqi_change_mean",
       sd_col = "abs_dlqi_change_sd",
       baseline_mean_col = "abs_dlqi_mean", baseline_sd_col = "abs_dlqi_sd",
       followup_mean_col = "abs_dlqi_mean", followup_sd_col = "abs_dlqi_sd",
       target_wk = 24),
  list(group = "safety", outcome = "sae",                view = "v_safety",
       kind = "binary", k_col = "sae",                target_wk = 24),
  list(group = "safety", outcome = "disc_any",           view = "v_safety",
       kind = "binary", k_col = "disc_any",           target_wk = 24),
  list(group = "safety", outcome = "disc_ae",            view = "v_safety",
       kind = "binary", k_col = "disc_ae",            target_wk = 24),
  list(group = "safety", outcome = "serious_infection",  view = "v_safety",
       kind = "binary", k_col = "serious_infection",  target_wk = 24),
  list(group = "safety", outcome = "injection_site_rxn", view = "v_safety",
       kind = "binary", k_col = "injection_site_rxn", target_wk = 24),
  list(group = "safety", outcome = "malignancy",         view = "v_safety",
       kind = "binary", k_col = "malignancy",         target_wk = 24),
  list(group = "safety", outcome = "nmsc",               view = "v_safety",
       kind = "binary", k_col = "nmsc",               target_wk = 24),
  list(group = "safety", outcome = "malignancy_non_nmsc", view = "v_safety",
       kind = "binary", k_col = "malignancy_non_nmsc", target_wk = 24)
)

# ---------------------------------------------------------------------------
# Pull each view once. The endpoints catalogue may consult the same view
# multiple times (e.g. v_safety has 8 outcomes) so cache.
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
# (ref_id, arm_no) with a non-NA value in `value_col`. Returns df subset.
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

# Same idea but for continuous endpoints: an arm contributes if EITHER its
# direct change is reported OR baseline + follow-up are both present (we
# derive the change). We still keep the latest timepoint per arm.
pick_primary_timepoint_continuous <- function(df, ep, target_wk) {
  has_change <- !is.na(df[[ep$mean_col]])
  has_bf     <- !is.na(df[[ep$followup_mean_col]])  # baseline backfilled later
  ok <- (has_change | has_bf) &
        !is.na(df$timepoint) & df$timepoint > 0 & df$timepoint <= target_wk
  d <- df[ok, , drop = FALSE]
  if (!nrow(d)) return(d[0, , drop = FALSE])
  key <- paste(d$ref_id, d$arm_no, sep = "|")
  ord <- order(key, -d$timepoint)
  d   <- d[ord, , drop = FALSE]
  d[!duplicated(paste(d$ref_id, d$arm_no, sep = "|")), , drop = FALSE]
}

# Mirror the app's derive_change(): when an arm reports baseline + follow-up
# but no direct change, fill in change_mean = follow - baseline; leave the
# SD missing (depends on unknown within-arm correlation).
attach_derived_change <- function(d, raw, ep) {
  if (!nrow(d)) {
    d$change_mean <- numeric(0); d$change_sd <- numeric(0); return(d)
  }
  # Baseline lookup from the raw view (timepoint == 0).
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

# ---------------------------------------------------------------------------
# Build the per-endpoint "long" table of contributing arms: one row per
# (ref_id, trial, arm_no, drug, n, k_or_mean[, sd]). The pairwise / single-
# arm meta-analyses all start from this shape.
# ---------------------------------------------------------------------------
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
# Per-comparison meta-analysis. For each (drug_a, drug_b) co-occurring in a
# trial, pull every trial's a-arm + b-arm and call metabin/metacont. If a
# trial has multiple arms of the same drug (e.g. dose variants) we sum n and
# k for binary outcomes, and pool the mean/sd via inverse-variance weighted
# average for continuous outcomes (rare; logged when it happens).
# ---------------------------------------------------------------------------
aggregate_arms_binary <- function(rows) {
  # rows: same trial, same drug, possibly multiple arms (dose variants).
  if (nrow(rows) == 1L) {
    return(list(n = rows$n[1], k = rows$k[1]))
  }
  list(n = sum(rows$n), k = sum(rows$k))
}

aggregate_arms_continuous <- function(rows) {
  if (nrow(rows) == 1L) {
    return(list(n = rows$n[1], mean = rows$mean[1], sd = rows$sd[1]))
  }
  # Inverse-variance weighted mean; pooled SD via the standard formula for
  # combining independent groups (treats the variants as a single arm).
  w  <- rows$n / (rows$sd^2)
  mu <- sum(w * rows$mean) / sum(w)
  # Pool variances around the combined mean.
  ss <- sum((rows$n - 1) * rows$sd^2 + rows$n * (rows$mean - mu)^2)
  v  <- ss / (sum(rows$n) - 1)
  list(n = sum(rows$n), mean = mu, sd = sqrt(v))
}

run_pairwise <- function(arm_tbl, drug_a, drug_b, ep) {
  ta <- arm_tbl[arm_tbl$drug == drug_a, , drop = FALSE]
  tb <- arm_tbl[arm_tbl$drug == drug_b, , drop = FALSE]
  trials_a <- unique(ta$ref_id)
  trials_b <- unique(tb$ref_id)
  trials   <- intersect(trials_a, trials_b)
  if (!length(trials)) return(NULL)

  per_trial <- vector("list", length(trials))
  for (i in seq_along(trials)) {
    rid <- trials[i]
    ra  <- ta[ta$ref_id == rid, , drop = FALSE]
    rb  <- tb[tb$ref_id == rid, , drop = FALSE]
    trial_name <- ra$trial[1] %||% rb$trial[1]
    if (ep$kind == "binary") {
      a <- aggregate_arms_binary(ra)
      b <- aggregate_arms_binary(rb)
      per_trial[[i]] <- data.frame(
        ref_id = rid, trial = trial_name,
        event_a = a$k, n_a = a$n,
        event_b = b$k, n_b = b$n,
        stringsAsFactors = FALSE
      )
    } else {
      a <- aggregate_arms_continuous(ra)
      b <- aggregate_arms_continuous(rb)
      per_trial[[i]] <- data.frame(
        ref_id = rid, trial = trial_name,
        mean_a = a$mean, sd_a = a$sd, n_a = a$n,
        mean_b = b$mean, sd_b = b$sd, n_b = b$n,
        stringsAsFactors = FALSE
      )
    }
  }
  pt <- do.call(rbind, per_trial)
  if (!nrow(pt)) return(NULL)

  ma <- tryCatch(suppressWarnings({
    if (ep$kind == "binary") {
      metabin(event.e = pt$event_a, n.e = pt$n_a,
              event.c = pt$event_b, n.c = pt$n_b,
              studlab = pt$trial,
              sm = "RR", method = "Inverse",
              method.tau = "REML",
              common = TRUE, random = TRUE,
              warn = FALSE)
    } else {
      metacont(n.e = pt$n_a, mean.e = pt$mean_a, sd.e = pt$sd_a,
               n.c = pt$n_b, mean.c = pt$mean_b, sd.c = pt$sd_b,
               studlab = pt$trial,
               sm = "MD", method.tau = "REML",
               common = TRUE, random = TRUE,
               warn = FALSE)
    }
  }), error = function(e) NULL)
  if (is.null(ma)) return(NULL)

  list(ma = ma, per_trial = pt)
}

# Expand meta::* object into the row we store.
ma_row <- function(ma, drug_a, drug_b, ep) {
  data.frame(
    drug_a         = drug_a,
    drug_b         = drug_b,
    endpoint_group = ep$group,
    outcome_code   = ep$outcome,
    sm             = ma$sm,
    n_studies      = ma$k,
    te_fe          = ma$TE.common,
    se_fe          = ma$seTE.common,
    lo_fe          = ma$lower.common,
    hi_fe          = ma$upper.common,
    z_fe           = ma$zval.common,
    p_fe           = ma$pval.common,
    te_re          = ma$TE.random,
    se_re          = ma$seTE.random,
    lo_re          = ma$lower.random,
    hi_re          = ma$upper.random,
    z_re           = ma$zval.random,
    p_re           = ma$pval.random,
    tau2           = ma$tau2 %||% NA_real_,
    i2             = ma$I2,
    q              = ma$Q,
    q_df           = ma$df.Q,
    q_pval         = ma$pval.Q,
    method_tau     = "REML",
    stringsAsFactors = FALSE
  )
}

trial_rows <- function(ma, pt, comparison_id, ep) {
  # meta stores per-study TE in $TE and SE in $seTE for both metabin and
  # metacont; weights in $w.common / $w.random.
  w_fe <- ma$w.common / sum(ma$w.common, na.rm = TRUE)
  w_re_sum <- sum(ma$w.random, na.rm = TRUE)
  w_re <- if (w_re_sum > 0) ma$w.random / w_re_sum else rep(NA_real_, length(ma$w.random))
  base <- data.frame(
    comparison_id = comparison_id,
    ref_id        = pt$ref_id,
    trial         = pt$trial,
    te            = ma$TE,
    se            = ma$seTE,
    lo            = ma$lower,
    hi            = ma$upper,
    weight_fe     = w_fe,
    weight_re     = w_re,
    stringsAsFactors = FALSE
  )
  if (ep$kind == "binary") {
    base$event_a <- pt$event_a; base$n_a <- pt$n_a
    base$event_b <- pt$event_b; base$n_b <- pt$n_b
    base$mean_a <- NA_real_; base$sd_a <- NA_real_
    base$mean_b <- NA_real_; base$sd_b <- NA_real_
  } else {
    base$event_a <- NA_integer_; base$n_a <- pt$n_a
    base$event_b <- NA_integer_; base$n_b <- pt$n_b
    base$mean_a <- pt$mean_a; base$sd_a <- pt$sd_a
    base$mean_b <- pt$mean_b; base$sd_b <- pt$sd_b
  }
  base
}

# ---------------------------------------------------------------------------
# Single-arm proportion meta-analysis (per drug, per binary outcome). Used
# by the node-filter view of the app.
# ---------------------------------------------------------------------------
run_proportion <- function(arm_tbl, drug, ep) {
  d <- arm_tbl[arm_tbl$drug == drug & !is.na(arm_tbl$k), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  # Collapse multi-arm (same drug, same trial) into a single entry. Using
  # aggregate with `by =` (rather than the formula form) so $k / $n stay
  # as plain numeric columns, not a matrix column.
  agg <- aggregate(d[, c("k", "n")],
                   by = list(ref_id = d$ref_id, trial = d$trial),
                   FUN = sum)
  if (!nrow(agg) || length(unique(agg$ref_id)) < 1) return(NULL)
  ma <- tryCatch(
    suppressWarnings(
      # method = "Inverse" forces classical inverse-variance pooling on the
      # logit scale; default ("GLMM") would otherwise force method.tau="ML".
      metaprop(event = agg$k, n = agg$n, studlab = agg$trial,
               sm = "PLOGIT", method = "Inverse", method.tau = "REML",
               common = TRUE, random = TRUE, warn = FALSE)
    ),
    error = function(e) NULL
  )
  if (is.null(ma)) return(NULL)
  list(ma = ma, per_trial = agg)
}

prop_row <- function(ma, drug, ep) {
  # Back-transform pooled logit to proportion + CI.
  inv <- function(x) exp(x) / (1 + exp(x))
  data.frame(
    drug           = drug,
    endpoint_group = ep$group,
    outcome_code   = ep$outcome,
    sm             = "proportion",
    n_studies      = ma$k,
    te_fe          = inv(ma$TE.common),
    lo_fe          = inv(ma$lower.common),
    hi_fe          = inv(ma$upper.common),
    te_re          = inv(ma$TE.random),
    lo_re          = inv(ma$lower.random),
    hi_re          = inv(ma$upper.random),
    tau2           = ma$tau2 %||% NA_real_,
    i2             = ma$I2,
    q              = ma$Q,
    q_df           = ma$df.Q,
    q_pval         = ma$pval.Q,
    method_tau     = "REML",
    stringsAsFactors = FALSE
  )
}

prop_trial_rows <- function(ma, agg, proportion_id) {
  w_fe <- ma$w.common / sum(ma$w.common, na.rm = TRUE)
  w_re_sum <- sum(ma$w.random, na.rm = TRUE)
  w_re <- if (w_re_sum > 0) ma$w.random / w_re_sum else rep(NA_real_, length(ma$w.random))
  # Wilson CI for the per-arm proportion (matches what metaprop reports).
  p_hat <- agg$k / agg$n
  z <- 1.959964
  denom <- 1 + z^2 / agg$n
  centre <- (p_hat + z^2 / (2 * agg$n)) / denom
  half   <- z * sqrt(p_hat * (1 - p_hat) / agg$n + z^2 / (4 * agg$n^2)) / denom
  data.frame(
    proportion_id = proportion_id,
    ref_id        = agg$ref_id,
    trial         = agg$trial,
    k             = agg$k,
    n             = agg$n,
    p             = p_hat,
    lo            = pmax(0, centre - half),
    hi            = pmin(1, centre + half),
    weight_fe     = w_fe,
    weight_re     = w_re,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Network meta-analysis. One NMA per (endpoint, outcome), using all drugs in
# the network. Returns:
#   - one summary row per network (k_studies, heterogeneity, inconsistency, ...)
#   - one row per ordered drug pair (drug_a vs drug_b) with FE and RE pooled
#     estimates, plus n_direct / n_indirect contributions.
#
# Multi-arm trials with the same drug at different doses are aggregated to a
# single arm per (trial, drug), matching the existing NMA connectivity diagram
# convention (nodes are drugs regardless of dose).
# ---------------------------------------------------------------------------
aggregate_per_trial_drug <- function(arm_tbl, ep) {
  key <- paste(arm_tbl$ref_id, arm_tbl$drug, sep = "|")
  if (!anyDuplicated(key)) return(arm_tbl)
  out <- split(arm_tbl, key, drop = TRUE)
  collapsed <- lapply(out, function(rows) {
    if (nrow(rows) == 1L) return(rows)
    if (ep$kind == "binary") {
      agg <- aggregate_arms_binary(rows)
      rows[1, "n"] <- agg$n
      rows[1, "k"] <- agg$k
    } else {
      agg <- aggregate_arms_continuous(rows)
      rows[1, "n"]    <- agg$n
      rows[1, "mean"] <- agg$mean
      rows[1, "sd"]   <- agg$sd
    }
    rows[1, , drop = FALSE]
  })
  do.call(rbind, collapsed)
}

# Build the long-format pairwise input that netmeta() expects. One row per
# direct comparison per trial. netmeta::pairwise() handles continuity
# correction for binary outcomes with zero events automatically.
build_nma_pairwise <- function(arm_tbl, ep) {
  arm_tbl <- aggregate_per_trial_drug(arm_tbl, ep)
  # Need at least 2 distinct treatments per study to form a comparison.
  by_trial <- split(arm_tbl, arm_tbl$ref_id)
  ok_ids <- names(by_trial)[vapply(by_trial,
                                   function(d) length(unique(d$drug)) >= 2L,
                                   logical(1))]
  if (!length(ok_ids)) return(NULL)
  arm_tbl <- arm_tbl[as.character(arm_tbl$ref_id) %in% ok_ids, , drop = FALSE]
  if (!nrow(arm_tbl)) return(NULL)
  args <- list(treat = arm_tbl$drug, studlab = arm_tbl$ref_id,
               data = arm_tbl, allstudies = TRUE)
  if (ep$kind == "binary") {
    args <- c(args, list(event = arm_tbl$k, n = arm_tbl$n, sm = "RR"))
  } else {
    args <- c(args, list(n = arm_tbl$n, mean = arm_tbl$mean, sd = arm_tbl$sd,
                         sm = "MD"))
  }
  pw <- tryCatch(suppressWarnings(do.call(meta::pairwise, args)),
                 error = function(e) NULL)
  if (is.null(pw) || !nrow(pw)) return(NULL)
  pw <- pw[is.finite(pw$TE) & is.finite(pw$seTE) & pw$seTE > 0, , drop = FALSE]
  if (!nrow(pw)) return(NULL)
  pw
}

# Run netmeta() on the pairwise long-format data. Returns the fitted object
# or NULL if the network is too small / disconnected.
run_nma <- function(pw) {
  treats <- unique(c(pw$treat1, pw$treat2))
  if (length(treats) < 2L) return(NULL)
  nm <- tryCatch(suppressWarnings(netmeta::netmeta(
    TE = pw$TE, seTE = pw$seTE,
    treat1 = pw$treat1, treat2 = pw$treat2,
    studlab = pw$studlab,
    data = pw,
    sm = pw$sm[1] %||% attr(pw, "sm") %||% "RR",
    common = TRUE, random = TRUE,
    reference.group = if ("Placebo" %in% treats) "Placebo" else treats[1],
    tol.multiarm = 0.5
  )), error = function(e) NULL)
  if (is.null(nm)) return(NULL)
  # Insist on a connected network.
  if (!isTRUE(nm$n.subnets == 1L) && !is.null(nm$n.subnets)) return(NULL)
  nm
}

# Extract one row per drug pair (drug_a vs drug_b). netmeta stores pooled
# TE / seTE / CI as treatment × treatment matrices; slice them into long form.
nma_pair_rows <- function(nm) {
  treats <- nm$trts
  n_t <- length(treats)
  pairs <- expand.grid(a = seq_len(n_t), b = seq_len(n_t),
                       KEEP.OUT.ATTRS = FALSE)
  pairs <- pairs[pairs$a != pairs$b, , drop = FALSE]
  # Direct evidence: pairwise long-format input has rows per direct study.
  pw <- nm$data
  count_direct <- function(ta, tb) {
    sum((pw$.treat1 == ta & pw$.treat2 == tb) |
        (pw$.treat1 == tb & pw$.treat2 == ta))
  }
  out <- data.frame(
    drug_a = treats[pairs$a],
    drug_b = treats[pairs$b],
    te_fe  = nm$TE.common[cbind(pairs$a, pairs$b)],
    se_fe  = nm$seTE.common[cbind(pairs$a, pairs$b)],
    lo_fe  = nm$lower.common[cbind(pairs$a, pairs$b)],
    hi_fe  = nm$upper.common[cbind(pairs$a, pairs$b)],
    te_re  = nm$TE.random[cbind(pairs$a, pairs$b)],
    se_re  = nm$seTE.random[cbind(pairs$a, pairs$b)],
    lo_re  = nm$lower.random[cbind(pairs$a, pairs$b)],
    hi_re  = nm$upper.random[cbind(pairs$a, pairs$b)],
    stringsAsFactors = FALSE
  )
  out$n_direct <- mapply(count_direct, out$drug_a, out$drug_b)
  # n_indirect: contributions to this pair from non-direct paths. netmeta
  # exposes contribution matrices but the simple n_studies - n_direct view is
  # a usable proxy; we expose total n_studies separately at the network level.
  out$n_indirect <- pmax(0L, nm$k - out$n_direct)
  out
}

# One-row network-level summary. q_total is decomposed into q_het (within-
# design) and q_inc (between-design / inconsistency).
nma_summary <- function(nm, ep, status = "ok") {
  data.frame(
    endpoint_group = ep$group,
    outcome_code   = ep$outcome,
    sm             = nm$sm,
    status         = status,
    n_studies      = nm$k,
    n_treatments   = length(nm$trts),
    n_pairwise     = nm$m,
    tau2           = nm$tau2 %||% NA_real_,
    i2             = nm$I2 %||% NA_real_,
    q_total        = nm$Q %||% NA_real_,
    q_het          = nm$Q.heterogeneity %||% NA_real_,
    q_inc          = nm$Q.inconsistency %||% NA_real_,
    p_inc          = nm$pval.Q.inconsistency %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

# Skeleton row when the network can't be fit (kept so the app can show a
# "too sparse" note rather than a missing block).
nma_sparse_summary <- function(ep) {
  data.frame(
    endpoint_group = ep$group,
    outcome_code   = ep$outcome,
    sm             = if (ep$kind == "binary") "RR" else "MD",
    status         = "sparse",
    n_studies      = NA_integer_,
    n_treatments   = NA_integer_,
    n_pairwise     = NA_integer_,
    tau2           = NA_real_,
    i2             = NA_real_,
    q_total        = NA_real_,
    q_het          = NA_real_,
    q_inc          = NA_real_,
    p_inc          = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Drive the full build.
# ---------------------------------------------------------------------------
all_pairwise        <- list()
all_pairwise_trials <- list()
all_proportion      <- list()
all_proportion_trials <- list()
all_nma             <- list()
all_nma_estimates   <- list()

cat("Building meta-analysis tables...\n")
for (ep in endpoints) {
  arm_tbl <- build_arm_table(ep)
  if (!nrow(arm_tbl)) {
    cat(sprintf("  %s / %s: no contributing arms - skipped\n",
                ep$group, ep$outcome))
    next
  }
  drugs <- sort(unique(arm_tbl$drug))

  # Pairwise: every co-occurring pair, and every (drug, Placebo) pair (the
  # latter may already be in the former; uniq removes duplicates).
  drug_pairs <- list()
  by_trial <- split(arm_tbl$drug, arm_tbl$ref_id)
  for (ds in by_trial) {
    ds <- sort(unique(ds))
    if (length(ds) < 2) next
    cmb <- utils::combn(ds, 2)
    for (j in seq_len(ncol(cmb))) {
      drug_pairs[[paste(cmb[1, j], cmb[2, j], sep = "|")]] <-
        c(cmb[1, j], cmb[2, j])
    }
  }
  # Force drug_a = "Placebo" goes to drug_b position so RR is "drug vs
  # placebo" (effect > 1 means drug helps for response outcomes).
  drug_pairs <- lapply(drug_pairs, function(p) {
    if (identical(p[1], "Placebo")) c(p[2], p[1]) else p
  })

  for (pair in drug_pairs) {
    res <- run_pairwise(arm_tbl, pair[1], pair[2], ep)
    if (is.null(res)) next
    row <- ma_row(res$ma, pair[1], pair[2], ep)
    cid <- sprintf("%s|%s|%s|%s",
                   ep$group, ep$outcome, pair[1], pair[2])
    row$comparison_id <- cid
    all_pairwise[[length(all_pairwise) + 1L]] <- row
    all_pairwise_trials[[length(all_pairwise_trials) + 1L]] <-
      trial_rows(res$ma, res$per_trial, cid, ep)
  }

  # Single-arm proportion (binary only).
  if (ep$kind == "binary") {
    for (drug in drugs) {
      pres <- run_proportion(arm_tbl, drug, ep)
      if (is.null(pres)) next
      prow <- prop_row(pres$ma, drug, ep)
      pid  <- sprintf("%s|%s|%s", ep$group, ep$outcome, drug)
      prow$proportion_id <- pid
      all_proportion[[length(all_proportion) + 1L]] <- prow
      all_proportion_trials[[length(all_proportion_trials) + 1L]] <-
        prop_trial_rows(pres$ma, pres$per_trial, pid)
    }
  }

  # Network meta-analysis for this endpoint × outcome. Same arm table as the
  # pairwise loop; doses collapsed to drug-level nodes.
  nma_status <- "skipped"
  pw <- build_nma_pairwise(arm_tbl, ep)
  if (is.null(pw)) {
    all_nma[[length(all_nma) + 1L]] <- nma_sparse_summary(ep)
    nma_status <- "sparse"
  } else {
    nm <- run_nma(pw)
    if (is.null(nm)) {
      all_nma[[length(all_nma) + 1L]] <- nma_sparse_summary(ep)
      nma_status <- "sparse"
    } else {
      ns <- nma_summary(nm, ep)
      nid <- sprintf("%s|%s", ep$group, ep$outcome)
      ns$network_id <- nid
      all_nma[[length(all_nma) + 1L]] <- ns
      ne <- nma_pair_rows(nm)
      ne$network_id <- nid
      all_nma_estimates[[length(all_nma_estimates) + 1L]] <- ne
      nma_status <- sprintf("ok (%d treats, %d studies)",
                            length(nm$trts), nm$k)
    }
  }
  # Tag the sparse skeleton row with a stable network_id too.
  if (identical(nma_status, "sparse")) {
    all_nma[[length(all_nma)]]$network_id <- sprintf("%s|%s",
                                                     ep$group, ep$outcome)
  }

  cat(sprintf("  %s / %s: %d pairwise, %d arms; NMA: %s\n",
              ep$group, ep$outcome,
              length(drug_pairs), length(unique(arm_tbl$drug)),
              nma_status))
}

bind_or_empty <- function(lst, template) {
  if (!length(lst)) template else do.call(rbind, lst)
}

ma_pairwise <- bind_or_empty(all_pairwise,
  data.frame(drug_a = character(), drug_b = character(),
             endpoint_group = character(), outcome_code = character(),
             stringsAsFactors = FALSE))
ma_pairwise_trials <- bind_or_empty(all_pairwise_trials,
  data.frame(comparison_id = character(), ref_id = integer(),
             trial = character(), stringsAsFactors = FALSE))
ma_proportion <- bind_or_empty(all_proportion,
  data.frame(drug = character(), endpoint_group = character(),
             outcome_code = character(), stringsAsFactors = FALSE))
ma_proportion_trials <- bind_or_empty(all_proportion_trials,
  data.frame(proportion_id = character(), ref_id = integer(),
             trial = character(), stringsAsFactors = FALSE))
ma_nma <- bind_or_empty(all_nma,
  data.frame(network_id = character(), endpoint_group = character(),
             outcome_code = character(), status = character(),
             stringsAsFactors = FALSE))
ma_nma_estimates <- bind_or_empty(all_nma_estimates,
  data.frame(network_id = character(), drug_a = character(),
             drug_b = character(), stringsAsFactors = FALSE))

built_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")
add_built_at <- function(df) {
  df$built_at <- if (nrow(df)) built_at else character(0)
  df
}
ma_pairwise          <- add_built_at(ma_pairwise)
ma_pairwise_trials   <- add_built_at(ma_pairwise_trials)
ma_proportion        <- add_built_at(ma_proportion)
ma_proportion_trials <- add_built_at(ma_proportion_trials)
ma_nma               <- add_built_at(ma_nma)
ma_nma_estimates     <- add_built_at(ma_nma_estimates)

cat(sprintf("\nWriting tables to %s ...\n", DB_PATH))
dbWriteTable(con, "ma_pairwise",          ma_pairwise,          overwrite = TRUE)
dbWriteTable(con, "ma_pairwise_trials",   ma_pairwise_trials,   overwrite = TRUE)
dbWriteTable(con, "ma_proportion",        ma_proportion,        overwrite = TRUE)
dbWriteTable(con, "ma_proportion_trials", ma_proportion_trials, overwrite = TRUE)
dbWriteTable(con, "ma_nma",               ma_nma,               overwrite = TRUE)
dbWriteTable(con, "ma_nma_estimates",     ma_nma_estimates,     overwrite = TRUE)

dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ma_pairwise_lookup
                ON ma_pairwise(endpoint_group, outcome_code, drug_a, drug_b)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ma_pairwise_trials_cid
                ON ma_pairwise_trials(comparison_id)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ma_proportion_lookup
                ON ma_proportion(endpoint_group, outcome_code, drug)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ma_proportion_trials_pid
                ON ma_proportion_trials(proportion_id)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ma_nma_lookup
                ON ma_nma(endpoint_group, outcome_code)")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ma_nma_estimates_lookup
                ON ma_nma_estimates(network_id, drug_a, drug_b)")

cat(sprintf("Done. ma_pairwise: %d, ma_pairwise_trials: %d, ma_proportion: %d, ma_proportion_trials: %d, ma_nma: %d, ma_nma_estimates: %d\n",
            nrow(ma_pairwise), nrow(ma_pairwise_trials),
            nrow(ma_proportion), nrow(ma_proportion_trials),
            nrow(ma_nma), nrow(ma_nma_estimates)))
