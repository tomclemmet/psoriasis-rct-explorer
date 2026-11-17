rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(multinma)
library(meta)
library(metafor)
options(mc.cores = parallel::detectCores())
load("C:/Users/p12916tc/Documents/Data/NMA/saved-results/meta-analysis.RData")
lookup <- read.csv("R/trt_class.csv")
con <- dbConnect(RSQLite::SQLite(), "pso-trial-explorer/psoriasis-rcts.sqlite")

# Extract data =================================================================

dbListTables(con)

pasi <- dbReadTable(con, "v_pasi")
dlqi <- dbReadTable(con, "v_dlqi")
safety <- dbReadTable(con, "v_safety")
join_keys <- colnames(pasi)[seq(1,9)]
dbDisconnect(con)

data <- pasi |> 
  full_join(dlqi, by = join_keys) |> 
  full_join(safety, by = join_keys) |> 
  filter(!is.na(drug)) |> 
  left_join(lookup, by = "drug")
  
drugs <- unique(data$drug)
comparisons <- as.data.frame(t(combn(drugs, 2)))
# results <- list()
niter <- 2000

# Network meta-analyses ========================================================
source("R/utils/app-process.R")
source("R/utils/load-pasi.R")
source("R/utils/ord-jags-nma.R")

## PASI Response ---------------------------------------------------------------

# JAGS .........................................................................

pasi <- readRDS("C:/Users/p12916tc/Documents/Data/NMA/saved-results/htai-pasi.rds")
results$pasi_fe_fez_ind_unadj <- nma_results(pasi$fe_fez_ind_unadj, label = "pasi", effects = "fixed", method = "fixed cutpoints")
results$pasi_re_fez_ind_unadj <- nma_results(pasi$re_fez_ind_unadj, label = "pasi", effects = "random", method = "fixed cutpoints")
results$pasi_fe_rezt_ind_unadj <- nma_results(pasi$fe_rezt_ind_unadj, label = "pasi", effects = "fixed", method = "treatment-level cutpoints")
results$pasi_re_rezt_ind_unadj <- nma_results(pasi$re_rezt_ind_unadj, label = "pasi", effects = "random", method = "treatment-level cutpoints")
results$pasi_fe_rezi_ind_unadj <- nma_results(pasi$fe_rezi_ind_unadj, label = "pasi", effects = "fixed", method = "trial-level cutpoints")
results$pasi_re_rezi_ind_unadj <- nma_results(pasi$re_rezi_ind_unadj, label = "pasi", effects = "random", method = "trial-level cutpoints")
results$pasi_fe_fez_ind_mu <- nma_results(pasi$fe_fez_ind_mu, label = "pasi", effects = "fixed", method = "fixed cutpoints & placebo adjustment")
results$pasi_re_fez_ind_mu <- nma_results(pasi$re_fez_ind_mu, label = "pasi", effects = "random", method = "fixed cutpoints & placebo adjustment")
results$pasi_fe_rezt_ind_mu <- nma_results(pasi$fe_rezt_ind_mu, label = "pasi", effects = "fixed", method = "treatment-level cutpoints & placebo adjustment")
results$pasi_re_rezt_ind_mu <- nma_results(pasi$re_rezt_ind_mu, label = "pasi", effects = "random", method = "treatment-level cutpoints & placebo adjustment")
results$pasi_fe_rezi_ind_mu <- nma_results(pasi$fe_rezi_ind_mu, label = "pasi", effects = "fixed", method = "trial-level cutpoints & placebo adjustment")
results$pasi_re_rezi_ind_mu <- nma_results(pasi$re_rezi_ind_mu, label = "pasi", effects = "random", method = "trial-level cutpoints & placebo adjustment")

dlqi <- readRDS("C:/Users/p12916tc/Documents/Data/NMA/saved-results/htai-dlqi.rds")
results$dlqi_fe_fez_ind_unadj <- nma_results(dlqi$fe_fez_ind_unadj, label = "dlqi", effects = "fixed", method = "fixed cutpoints")
results$dlqi_re_fez_ind_unadj <- nma_results(dlqi$re_fez_ind_unadj, label = "dlqi", effects = "random", method = "fixed cutpoints")
results$dlqi_fe_rezt_ind_unadj <- nma_results(dlqi$fe_rezt_ind_unadj, label = "dlqi", effects = "fixed", method = "treatment-level cutpoints")
results$dlqi_re_rezt_ind_unadj <- nma_results(dlqi$re_rezt_ind_unadj, label = "dlqi", effects = "random", method = "treatment-level cutpoints")
results$dlqi_fe_rezi_ind_unadj <- nma_results(dlqi$fe_rezi_ind_unadj, label = "dlqi", effects = "fixed", method = "trial-level cutpoints")
results$dlqi_re_rezi_ind_unadj <- nma_results(dlqi$re_rezi_ind_unadj, label = "dlqi", effects = "random", method = "trial-level cutpoints")
results$dlqi_fe_fez_ind_mu <- nma_results(dlqi$fe_fez_ind_mu, label = "dlqi", effects = "fixed", method = "fixed cutpoints & placebo adjustment")
results$dlqi_re_fez_ind_mu <- nma_results(dlqi$re_fez_ind_mu, label = "dlqi", effects = "random", method = "fixed cutpoints & placebo adjustment")
results$dlqi_fe_rezt_ind_mu <- nma_results(dlqi$fe_rezt_ind_mu, label = "dlqi", effects = "fixed", method = "treatment-level cutpoints & placebo adjustment")
results$dlqi_re_rezt_ind_mu <- nma_results(dlqi$re_rezt_ind_mu, label = "dlqi", effects = "random", method = "treatment-level cutpoints & placebo adjustment")
results$dlqi_fe_rezi_ind_mu <- nma_results(dlqi$fe_rezi_ind_mu, label = "dlqi", effects = "fixed", method = "trial-level cutpoints & placebo adjustment")
results$dlqi_re_rezi_ind_mu <- nma_results(dlqi$re_rezi_ind_mu, label = "dlqi", effects = "random", method = "trial-level cutpoints & placebo adjustment")


## Absolute change in PASI -----------------------------------------------------

abs_pasi_data <- data |> 
  filter(if_any(contains("abs_pasi"), \(x) !is.na(x))) |> 
  select(ref_id, arm_no, n, drug, timepoint, contains("abs_pasi")) |> 
  mutate(baseline = if_else(timepoint == 0, "baseline", "follow_up")) |> 
  pivot_wider(names_from = baseline, 
              values_from = c(timepoint, abs_pasi_mean, abs_pasi_sd)) |> 
  mutate(
    abs_pasi_change_mean = if_else(
      is.na(abs_pasi_change_mean), 
      abs_pasi_mean_follow_up - abs_pasi_mean_baseline, abs_pasi_change_mean
    ),
    abs_pasi_change_sd = if_else(
      is.na(abs_pasi_change_sd), 
      sqrt((abs_pasi_sd_follow_up)^2 + (abs_pasi_sd_baseline)^2 - 
             2 * 0.5 * abs_pasi_sd_follow_up * abs_pasi_sd_baseline), 
      abs_pasi_change_sd
    )
  ) |> 
  filter(!is.na(abs_pasi_change_mean) & !is.na(abs_pasi_change_sd))

abs_pasi_ref <- metagen(
  TE = abs_pasi_change_mean, 
  seTE = abs_pasi_change_sd / sqrt(n), 
  data = filter(abs_pasi_data, drug == "Placebo")
)

abs_pasi_net <- set_agd_arm(
  abs_pasi_data, 
  study = ref_id,
  trt = drug,
  y = abs_pasi_change_mean, 
  se = abs_pasi_change_sd / sqrt(n),
  sample_size = n,
  trt_ref = "Placebo"
)

abs_pasi_fit_fe <- nma(
  abs_pasi_net,
  trt_effects = "fixed",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  iter = niter
)

results$abs_pasi_fe <- nma_results(
  abs_pasi_fit_fe, 
  distr(qnorm, abs_pasi_ref$TE.random, abs_pasi_ref$seTE.random),
  label = "abs_pasi_change"
)

abs_pasi_fit_re <- nma(
  abs_pasi_net,
  trt_effects = "random",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  iter = niter
)

results$abs_pasi_re <- nma_results(
  abs_pasi_fit_re, 
  distr(qnorm, abs_pasi_ref$TE.random, abs_pasi_ref$seTE.random),
  label = "abs_pasi_change"
)

## Absolute change in DLQI -----------------------------------------------------

abs_dlqi_data <- data |> 
  filter(if_any(contains("abs_dlqi"), \(x) !is.na(x))) |> 
  select(ref_id, arm_no, n, drug, timepoint, contains("abs_dlqi")) |> 
  mutate(baseline = if_else(timepoint == 0, "baseline", "follow_up")) |> 
  pivot_wider(names_from = baseline, 
              values_from = c(timepoint, abs_dlqi_mean, abs_dlqi_sd)) |> 
  mutate(
    abs_dlqi_change_mean = if_else(
      is.na(abs_dlqi_change_mean), 
      abs_dlqi_mean_follow_up - abs_dlqi_mean_baseline, abs_dlqi_change_mean
    ),
    abs_dlqi_change_sd = if_else(
      is.na(abs_dlqi_change_sd), 
      # Assumed 0.5 covariance
      sqrt((abs_dlqi_sd_follow_up)^2 + (abs_dlqi_sd_baseline)^2 - 
             2 * 0.5 * abs_dlqi_sd_follow_up * abs_dlqi_sd_baseline), 
      abs_dlqi_change_sd 
    )
  ) |> 
  filter(!is.na(abs_dlqi_change_mean) & !is.na(abs_dlqi_change_sd))

abs_dlqi_ref <- metagen(
  TE = abs_dlqi_change_mean, 
  seTE = abs_dlqi_change_sd / sqrt(n), 
  data = filter(abs_dlqi_data, drug == "Placebo")
)

abs_dlqi_net <- set_agd_arm(
  abs_dlqi_data, 
  study = ref_id,
  trt = drug,
  y = abs_dlqi_change_mean, 
  se = abs_dlqi_change_sd / sqrt(n),
  sample_size = n,
  trt_ref = "Placebo"
)

abs_dlqi_fit_fe <- nma(
  abs_dlqi_net,
  trt_effects = "fixed",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  iter = niter
)

results$abs_dlqi_fe <- nma_results(
  abs_dlqi_fit_fe, 
  distr(qnorm, abs_dlqi_ref$TE.random, abs_dlqi_ref$seTE.random),
  label = "abs_dlqi_change"
)

abs_dlqi_fit_re <- nma(
  abs_dlqi_net,
  trt_effects = "random",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  iter = niter
)

results$abs_dlqi_re <- nma_results(
  abs_dlqi_fit_re, 
  distr(qnorm, abs_dlqi_ref$TE.random, abs_dlqi_ref$seTE.random),
  label = "abs_dlqi_change"
)

## Binary outcomes -------------------------------------------------------------

bin_outcomes <- c(
  "pasi50", "pasi75", "pasi90", "pasi100", "dlqi_0_1", "dlqi_0", "sae", 
  "disc_any", "disc_ae", "serious_infection", "injection_site_rxn", "malignancy"
)
bin_outcomes <- c("dlqi_0_1", "dlqi_0")

bin_fit_fe <- list()
bin_fit_re <- list()

for (i in 1:length(bin_outcomes)) {
  
  df <- data |> 
    group_by(ref_id) |> 
    filter(all(!is.na(.data[[bin_outcomes[i]]]))) |> 
    ungroup()
  
  bin_net <- set_agd_arm(
    df,
    study = ref_id,
    trt = drug,
    r = .data[[bin_outcomes[i]]],
    n = n,
    trt_ref = "Placebo"
  )
  
  bin_fit_fe[[i]] <- nma(
    bin_net,
    trt_effects = "fixed",
    prior_intercept = normal(scale = 100),
    prior_trt = normal(scale = 10),
    iter = niter
  )
  
  bin_fit_re[[i]] <- nma(
    bin_net,
    trt_effects = "random",
    prior_intercept = normal(scale = 100),
    prior_trt = normal(scale = 10),
    prior_het = half_normal(scale = 5),
    iter = niter
  )
  
  message(bin_outcomes[i])
}

for (i in 1:length(bin_outcomes)) {
  placebo_data <- data |> 
    filter(drug == "Placebo", !is.na(.data[[bin_outcomes[i]]]))
    
  bin_ref <- metaprop(
    event = placebo_data[[bin_outcomes[i]]],
    n = n,
    data = placebo_data
  )
  
  results[[paste(bin_outcomes[i], "fe")]] <- nma_results(
    bin_fit_fe[[i]],
    beta_dist_metaprop(bin_ref, "random"),
    label = bin_outcomes[i]
  )

  results[[paste(bin_outcomes[i], "re")]] <- nma_results(
    bin_fit_re[[i]],
    beta_dist_metaprop(bin_ref, "random"),
    label = bin_outcomes[i]
  )

  message(bin_outcomes[i])
}

# Pairwise Meta-Analyses =======================================================

outcomes <- c(
  "pasi50", "pasi75", "pasi90", "pasi100",
  "dlqi_0_1", "dlqi_0", 
  "sae", "disc_any", "disc_ae", "serious_infection", "injection_site_rxn", 
  "malignancy"
)

## Binary outcomes -------------------------------------------------------------

for (i in 1:length(outcomes)) {
  for (j in 1:nrow(comparisons)) {
    tx <- comparisons[[j, 1]]
    ref <- comparisons[[j, 2]]
    outcome <- outcomes[i]
    comp_data <- data |> 
      group_by(ref_id) |> 
      filter(any(drug == tx) & any(drug == ref),
             drug %in% c(tx, ref)) |> 
      ungroup() |> 
      mutate(drug = if_else(drug == tx, "tx", "ref")) |> 
      select(ref_id, arm_no, drug, n, contains(outcome)) |> 
      filter(!is.na(.data[[outcome]])) |> 
      summarise(
        .by = c(ref_id, drug),
        k = sum(.data[[outcome]]),
        n = sum(n),
      ) |> 
      pivot_wider(names_from = drug, values_from = c(n, k))
    
    if(nrow(comp_data) <= 1) next
      
    fit <- metabin(
      event.e = comp_data$k_tx, n.e = comp_data$n_tx,event.c = comp_data$k_ref,
      n.c = comp_data$n_ref, sm = "RD"
    )
    results[[paste(outcome, tx, ref)]] <- nma_results(
      fit, label = outcome, t = tx, reft = ref
    )
  }
  message(outcome)
}

## Absolute change in PASI -----------------------------------------------------

for (j in 1:nrow(comparisons)) {
  tx <- comparisons[[j, 1]]
  ref <- comparisons[[j, 2]]
  pairwise <- abs_pasi_data |> 
    group_by(ref_id) |> 
    filter(any(drug == tx) & any(drug == ref),
           drug %in% c(tx, ref)) |> 
    ungroup() |> 
    mutate(drug = if_else(drug == tx, "tx", "ref")) |> 
    summarise(
      .by = c(ref_id, drug),
      n = sum(n),
      mu = mean(abs_pasi_change_mean),
      sd = mean(abs_pasi_change_sd)
    ) |> 
    pivot_wider(names_from = drug, values_from = c(mu, sd, n))
  
  if (nrow(pairwise) <= 1) next
  
  fit <- metacont(
    n.e = pairwise$n_tx, mean.e = pairwise$mu_tx, sd.e = pairwise$sd_tx,
    n.c = pairwise$n_ref, mean.c = pairwise$mu_ref, sd.c = pairwise$sd_ref,
    studlab = pairwise$ref_id, sm = "MD"
  )
  
  results[[paste("abs_pasi_change", tx, ref)]] <- nma_results(
    fit, label = "abs_pasi_change", t = tx, reft = ref
  )
}

## Absolute change in DLQI -----------------------------------------------------

for (j in 1:nrow(comparisons)) {
  tx <- comparisons[[j, 1]]
  ref <- comparisons[[j, 2]]
  pairwise <- abs_dlqi_data |> 
    group_by(ref_id) |> 
    filter(any(drug == tx) & any(drug == ref),
           drug %in% c(tx, ref)) |> 
    ungroup() |> 
    mutate(drug = if_else(drug == tx, "tx", "ref")) |> 
    summarise(
      .by = c(ref_id, drug),
      n = sum(n),
      mu = mean(abs_dlqi_change_mean),
      sd = mean(abs_dlqi_change_sd)
    ) |> 
    pivot_wider(names_from = drug, values_from = c(mu, sd, n))
  
  if (nrow(pairwise) <= 1) next
  
  fit <- metacont(
    n.e = pairwise$n_tx, mean.e = pairwise$mu_tx, sd.e = pairwise$sd_tx,
    n.c = pairwise$n_ref, mean.c = pairwise$mu_ref, sd.c = pairwise$sd_ref,
    studlab = pairwise$ref_id, sm = "MD"
  )
  
  results[[paste("abs_dlqi_change", tx, ref)]] <- nma_results(
    fit, label = "abs_dlqi_change", t = tx, reft = ref
  )
}

# Univariate meta-analysis =====================================================

drugs <- unique(data$drug)

## Binary outcomes -------------------------------------------------------------

for (i in 1:length(outcomes)) {
  for (k in 1:length(drugs)) {
    univar <- data |> 
      filter(
        drug == drugs[k],
        !is.na(.data[[outcomes[i]]])
      )
    
    if(nrow(univar) <= 1) next
    
    fit <- metaprop(
      univar[[outcomes[i]]], 
      univar$n, 
      studylab = univar$ref_id
    )
    results[[paste(outcomes[i], drugs[k])]] <- nma_results(
      fit, label = outcomes[i], t = drugs[k]
    )
  }
}

## Absolute change in PASI -----------------------------------------------------

for (k in 1:length(drugs)) {
  univar <- abs_pasi_data |> 
    filter(drug == drugs[k])
  
  if(nrow(univar) <= 1) next
  
  fit <- metagen(
    TE = univar$abs_pasi_change_mean, 
    seTE = univar$abs_pasi_change_sd / sqrt(univar$n),
    studylab = univar$ref_id
  )
  
  results[[paste("abs_pasi_change", drugs[k])]] <- nma_results(
    fit, label = "abs_pasi_change", t = drugs[k]
  )
}

## Absolute change in DLQI -----------------------------------------------------

for (k in 1:length(drugs)) {
  univar <- abs_dlqi_data |> 
    filter(drug == drugs[k])
  
  if(nrow(univar) <= 1) next
  
  fit <- metagen(
    TE = univar$abs_dlqi_change_mean, 
    seTE = univar$abs_dlqi_change_sd / sqrt(univar$n),
    studylab = univar$ref_id
  )
  
  results[[paste("abs_dlqi_change", drugs[k])]] <- nma_results(
    fit, label = "abs_dlqi_change", t = drugs[k]
  )
}

# Write results ================================================================
con <- dbConnect(RSQLite::SQLite(), "pso-trial-explorer/psoriasis-rcts.sqlite")

results_df <- bind_rows(results) |> 
  mutate(endpoint_group = if_else(likelihood == "multinomial" & 
                                    str_detect(endpoint, "pasi|dlqi"),
                                  substr(endpoint, 1, 4), endpoint)) |> 
  arrange(type != "network", likelihood != "multinomial") |> 
  group_by(endpoint_group, type, likelihood, method, effects, dic) |> 
  mutate(ma_id = cur_group_id()) |> ungroup()

model_info <- results_df |> 
  distinct(ma_id, .keep_all = TRUE) |> 
  mutate(metareg = if_else(likelihood == "multinomial", coalesce(str_extract(method, "(?<= & ).*$"), "unadjusted"), "unadjusted"),
         method = if_else(likelihood == "multinomial", coalesce(str_extract(method, "^.*(?= & )"), method), method))

model_estimates <- results_df |> 
  select(-c(type, likelihood, method, effects, dic, endpoint_group))

dbWriteTable(con, name = "ma_models", value = model_info, overwrite = TRUE)
dbWriteTable(con, name = "ma_results", value = model_estimates, overwrite = TRUE)

create_view_sql <- "
  CREATE VIEW v_meta_analysis AS
  SELECT *
  FROM ma_results
  LEFT JOIN ma_models ON ma_results.ma_id = ma_models.ma_id
"
dbExecute(con, "DROP VIEW IF EXISTS v_meta_analysis")
dbExecute(con, create_view_sql)

dbDisconnect(con)
save.image("C:/Users/p12916tc/Documents/Data/NMA/saved-results/meta-analysis.RData")
