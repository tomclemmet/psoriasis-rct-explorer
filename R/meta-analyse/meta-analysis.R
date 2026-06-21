rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)
library(meta)
options(mc.cores = parallel::detectCores())
source("R/meta-analyse/ma-utils.R")

# Extract data =================================================================

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")
dbListTables(con)

query <- "
SELECT 
    p.ref_id, p.arm_no, p.drug, p.timepoint, p.timepoint_unit, p.n,
    p.pasi50, p.pasi75, p.pasi90, p.pasi100, p.abs_pasi_change_mean, 
    p.abs_pasi_change_sd, p.abs_pasi_mean, p.abs_pasi_sd,
    d.dlqi_0_1, d.dlqi_0, d.abs_dlqi_change_mean, d.abs_dlqi_change_sd,
    d.abs_dlqi_mean, d.abs_dlqi_sd,
    s.sae, s.disc_any, s.disc_ae, s.serious_infection, s.injection_site_rxn,
    s.malignancy
FROM v_pasi p
LEFT JOIN v_dlqi d      ON p.ref_id = d.ref_id   AND p.arm_no = d.arm_no
                       AND p.timepoint = d.timepoint     
LEFT JOIN v_safety s    ON p.ref_id = s.ref_id   AND p.arm_no = s.arm_no
                       AND p.timepoint = s.timepoint
"

data <- dbGetQuery(con, query) |>
  filter(!is.na(drug)) |> 
  group_by(ref_id, arm_no) |> 
  filter(timepoint == 0 | timepoint == max(timepoint)) |> 
  ungroup()
drugs <- unique(data$drug)
comparisons <- as.data.frame(t(combn(drugs, 2)))
results <- list()
results <- readRDS("R/meta-analyse/ma-results.rds")
niter <- 2000

# Network meta-analyses ========================================================

## PASI Response ---------------------------------------------------------------

pasi_net <- set_agd_arm(
  filter(data, !if_all(pasi50:pasi100, \(x) is.na(x))),
  study = ref_id,
  trt = drug,
  r =  multi(r0 = n,
             pasi50, pasi75, pasi90, pasi100,
             inclusive = FALSE,
             type = "ordered")
)
pasi_fit_fe <- nma(
  pasi_net,
  trt_effects = "fixed",
  link = "probit",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  prior_aux = flat(),
  iter = niter
)

pasi_ref_fe <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(pasi_fit_fe, pars = "mu"))
)

results$pasi_fe <- nma_results(
  pasi_fit_fe, 
  pasi_ref_fe$TE.fixed, 
  pasi_ref_fe$seTE.fixed
)

# Random effects
pasi_fit_re <- nma(
  pasi_net,
  trt_effects = "random",
  link = "probit",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  prior_aux = flat(),
  iter = niter
)

pasi_ref_re <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(pasi_fit_re, pars = "mu"))
)

results$pasi_re <- nma_results(
  pasi_fit_re, 
  pasi_ref_re$TE.random, 
  pasi_ref_re$seTE.random
)

## DLQI response ---------------------------------------------------------------

dlqi_net <- set_agd_arm(
  filter(data, !if_all(dlqi_0_1:dlqi_0, \(x) is.na(x))),
  study = ref_id,
  trt = drug,
  r =  multi(r0 = n,
             dlqi_0_1, dlqi_0,
             inclusive = FALSE,
             type = "ordered")
)
dlqi_fit_fe <- nma(
  dlqi_net,
  trt_effects = "fixed",
  link = "probit",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  prior_aux = flat(),
  iter = niter
)

dlqi_ref_fe <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(dlqi_fit_fe, pars = "mu"))
)

results$dlqi_fe <- nma_results(
  dlqi_fit_fe, 
  dlqi_ref_fe$TE.fixed, 
  dlqi_ref_fe$seTE.fixed
)

dlqi_fit_re <- nma(
  dlqi_net,
  trt_effects = "random",
  link = "probit",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  prior_aux = flat(),
  iter = niter
)

dlqi_ref_re <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(dlqi_fit_re, pars = "mu"))# Is it correct to average mu?
)

results$dlqi_re <- nma_results(
  dlqi_fit_re, 
  dlqi_ref_re$TE.random, 
  dlqi_ref_re$seTE.random
)

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
      sqrt((abs_pasi_sd_follow_up)^2 + (abs_pasi_sd_baseline)^2 - 2 * 0.5 * abs_pasi_sd_follow_up * abs_pasi_sd_baseline), 
      abs_pasi_change_sd
    )
  ) |> 
  filter(!is.na(abs_pasi_change_mean) & !is.na(abs_pasi_change_sd))

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

abs_pasi_ref_fe <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(abs_pasi_fit_fe, pars = "mu"))# Is it correct to average mu?
)

results$abs_pasi_fe <- nma_results(
  abs_pasi_fit_fe, 
  abs_pasi_ref_fe$TE.fixed, 
  abs_pasi_ref_fe$seTE.fixed,
  label = "abs_pasi_change"
)

abs_pasi_fit_re <- nma(
  abs_pasi_net,
  trt_effects = "random",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  iter = niter
)

abs_pasi_ref_re <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(abs_pasi_fit_re, pars = "mu"))# Is it correct to average mu?
)

results$abs_pasi_re <- nma_results(
  abs_pasi_fit_re, 
  abs_pasi_ref_re$TE.random, 
  abs_pasi_ref_re$seTE.random,
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
      sqrt((abs_dlqi_sd_follow_up)^2 + (abs_dlqi_sd_baseline)^2 - 2 * 0.5 * abs_dlqi_sd_follow_up * abs_dlqi_sd_baseline), abs_dlqi_change_sd # Assumed 0.5 covariance
    )
  ) |> 
  filter(!is.na(abs_dlqi_change_mean) & !is.na(abs_dlqi_change_sd))

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

abs_dlqi_ref_fe <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(abs_dlqi_fit_fe, pars = "mu"))# Is it correct to average mu?
)

results$abs_dlqi_fe <- nma_results(
  abs_dlqi_fit_fe, 
  abs_dlqi_ref_fe$TE.fixed, 
  abs_dlqi_ref_fe$seTE.fixed,
  label = "abs_dlqi_change"
)

abs_dlqi_fit_re <- nma(
  abs_dlqi_net,
  trt_effects = "random",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  iter = niter
)

abs_dlqi_ref_re <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(abs_dlqi_fit_re, pars = "mu"))# Is it correct to average mu?
)

results$abs_dlqi_re <- nma_results(
  abs_dlqi_fit_re, 
  abs_dlqi_ref_re$TE.random, 
  abs_dlqi_ref_re$seTE.random,
  label = "abs_dlqi_change"
)

## Binary outcomes -------------------------------------------------------------

bin_outcomes <- c(
  "sae", "disc_any", "disc_ae", "serious_infection", "injection_site_rxn", 
  "malignancy"
)

for (i in 1:length(bin_outcomes)) {
  bin_net <- set_agd_arm(
    filter(data, !is.na(.data[[bin_outcomes[i]]])),
    study = ref_id,
    trt = drug,
    r = .data[[bin_outcomes[i]]],
    n = n,
    trt_ref = "Placebo"
  )
  
  bin_fit_fe <- nma(
    bin_net,
    trt_effects = "fixed",
    prior_intercept = normal(scale = 100),
    prior_trt = normal(scale = 10),
    iter = niter
  )
  
  bin_ref_fe <- metagen(
    TE = mean, 
    seTE = sd, 
    data = as.data.frame(summary(bin_fit_fe, pars = "mu"))
  )
  
  results[[paste(bin_outcomes[i], "fe")]] <- nma_results(
    bin_fit_fe, 
    bin_ref_fe$TE.fixed, 
    bin_ref_fe$seTE.fixed,
    label = bin_outcomes[i]
  )
  
  # Random effects
  bin_fit_re <- nma(
    bin_net,
    trt_effects = "random",
    prior_intercept = normal(scale = 100),
    prior_trt = normal(scale = 10),
    prior_het = half_normal(scale = 5),
    iter = niter
  )
  
  bin_ref_re <- metagen(
    TE = mean, 
    seTE = sd, 
    data = as.data.frame(summary(bin_fit_re, pars = "mu"))
  )
  
  results[[paste(bin_outcomes[i], "re")]] <- nma_results(
    bin_fit_re, 
    bin_ref_re$TE.random, 
    bin_ref_re$seTE.random,
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
    
    fit <- metaprop(univar[[outcomes[i]]], univar$n, studylab = univar$ref_id)
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
  
  results[[paste("abs_change_pasi", drugs[k])]] <- nma_results(
    fit, label = "abs_change_pasi", t = drugs[k]
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
  
  results[[paste("abs_change_dlqi", drugs[k])]] <- nma_results(
    fit, label = "abs_change_dlqi", t = drugs[k]
  )
}

# Write results ================================================================

results_table <- bind_rows(results)

dbWriteTable(con, name = "meta_analysis", value = results_table, overwrite = TRUE)

create_view_sql <- "
  CREATE VIEW v_meta_analysis AS
  SELECT *
  FROM meta_analysis
"
dbExecute(con, create_view_sql)

dbDisconnect(con)
saveRDS(results, "R/meta-analyse/ma-results.rds")
