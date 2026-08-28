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
load("R/meta-analyse/meta-analysis.RData")
lookup <- read.csv("R/meta-analyse/trt_class.csv")
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

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
results <- list()
niter <- 2000

pasi_drugs <- c("Placebo", setdiff(sort(unique(pasi$drug)), "Placebo"))

# Network meta-analyses ========================================================
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")

## PASI Response ---------------------------------------------------------------

# JAGS .........................................................................

pasi_ref <- metaprop(
  event = pasi50,
  n = n,
  data = filter(data, drug == "Placebo", !is.na(pasi50)),
  # sm = "PLOGIT",
  # method = "Inverse",
  # method.incr = "all",
  # incr = 0.5
)

settings <- expand.grid(
  effects = c("fixed", "random"), cutpoints = c("fixed", "random"), 
  baseline = c("unadjusted", "adjusted"), class = c("independent", "exchangeable"),
  consistency = c("consistency", "ume"), stringsAsFactors = FALSE) |> 
  filter(!(class == "exchangeable" & consistency == "ume")) |> 
  mutate(filename = paste(if_else(effects == "fixed", "fe", "re"),
                          if_else(cutpoints == "fixed", "fez", "rez"),
                          if_else(baseline == "unadjusted", "u", "a"),
                          if_else(class == "independent", "nc", "c"),
                          if_else(consistency == "consistency", "con", "ume"),
                          sep = "_"))

for (i in 1:nrow(settings)) {
  message(paste0("Model ", i, " of " , nrow(settings)))
  j[[settings$filename[i]]] <- pso_jags(
    pasi_jags, 
    filename = paste0("JAGS/", settings$filename[i], ".jags"), 
    effects = settings$effects[i], 
    cutpoints = settings$cutpoints[i], 
    baseline = settings$baseline[i],
    class = settings$class[i],
    consistency = settings$consistency[i]
  )
}

mod <- pso_jags(pasi_jags, consistency = "ume")
process_jags(j$re_rez_a_c_ume)

compare_jags(j)
devplot(j$re_rez_u_nc_con, j$re_rez_a_nc_con, output = "plot")
lapply(j, \(x) {process_jags(x)$DIC}) |> as.data.frame()

j |>
  lapply(\(x) {process_jags(x)$summary}) |>
  bind_rows(.id = "id") |>
  filter(!is.na(drug), ! drug %in% c("Phototherapy", "Mirikizumab", "Placebo",
                                     "Izokibep")) |>
  left_join(lookup, by = "drug") |>
  mutate(drug = forcats::fct_reorder(drug, mean, .fun = base::mean)) |>
  ggplot(aes(x = mean, y = drug, colour = id)) +
  geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`),
                  position = position_dodge(width = 0.7), shape = 15,
                  size = 0.1) +
  scale_colour_viridis_d(option = "turbo") +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_wrap(~ class, scales = "free_y")
ggsave("output/forest.png", height = 7, width = 10)

drug_rank <- process_jags(j$re_rez_a_c_con)$summary |> 
  filter(str_starts(param, "prob")) |> 
  mutate(drug = pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))]) |> 
  slice_head(n = 1, by = drug) |> 
  arrange(mean)

process_jags(j$re_rez_a_c_con)$summary |> 
  filter(str_starts(param, "prob")) |> 
  mutate(
    drug = factor(
      pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))],
      levels = drug_rank$drug
    ),
    outcome = factor(
      outcomes[as.numeric(str_extract(param, "(?<=\\[).*?(?=,)"))],
      levels = c("pasi50", "pasi75", "pasi90", "pasi100"),
      labels = c("PASI 50-75", "PASI 75-90", "PASI 90-100", "PASI 100")
    )
  ) |> 
  arrange(drug, desc(outcome)) |> 
  mutate(.by = drug, mean = mean - lag(mean, default = 0), .after = mean) |> 
  mutate(outcome = forcats::fct_rev(outcome)) |> 
  filter(drug %notin% exc) |> 
  ggplot(aes(x = mean, y = drug)) +
  geom_col(aes(fill = outcome), position = position_stack(reverse = TRUE)) +
  theme_classic() +
  scale_fill_viridis_d() +
  theme(legend.position = "bottom") +
  labs(title = "RE REZ with baseline adjustment")
ggsave("output/props.png", height = 7, width = 4)
  

results$fe_fez_u <- nma_results(j$fe_fez_u_nc_con, 
                                effects = "fixed", 
                                method = "standard")
results$re_fez_u <- nma_results(j$re_fez_u_nc_con, 
                                effects = "random", 
                                method = "standard")
results$fe_rez_u <- nma_results(j$fe_rez_u_nc_con, 
                                effects = "fixed", 
                                method = "REZ")
results$re_rez_u <- nma_results(j$re_rez_u_nc_con, 
                                effects = "random",
                                method = "REZ")
results$fe_fez_a <- nma_results(j$fe_fez_a_nc_con, 
                                effects = "fixed", 
                                method = "baseline adjusted")
results$re_fez_a <- nma_results(j$re_fez_a_nc_con, 
                                effects = "random", 
                                method = "baseline adjusted")
results$fe_rez_a <- nma_results(j$fe_rez_a_nc_con, 
                                effects = "fixed", 
                                method = "REZ, baseline adjusted")
results$re_rez_a <- nma_results(j$re_rez_a_nc_con, 
                                effects = "random", 
                                method = "REZ, baseline adjusted")

## multinma ....................................................................
# 
# pasi_net <- set_agd_arm(
#   filter(data, !if_all(pasi50:pasi100, \(x) is.na(x))),
#   study = ref_id,
#   trt = drug,
#   r =  multi(r0 = n,
#              pasi50, pasi75, pasi90, pasi100,
#              inclusive = TRUE,
#              type = "ordered")
#   )
# 
# pasi_fit_fe <- nma(
#   pasi_net,
#   trt_effects = "fixed",
#   link = "probit",
#   prior_intercept = normal(scale = 100),
#   prior_trt = normal(scale = 10),
#   prior_aux = flat(),
#   iter = 500,
#   chains = 2
# )
# 
# pasi_fit_fe_nodesplit <- nma(
#   pasi_net,
#   consistency = "nodesplit",
#   trt_effects = "fixed",
#   link = "probit",
#   prior_intercept = normal(scale = 100),
#   prior_trt = normal(scale = 10),
#   prior_aux = flat(),
#   iter = 500,
#   chains = 2
# )
# 
# pasi_fit_fe_ume <- nma(
#   pasi_net,
#   consistency = "ume",
#   trt_effects = "fixed",
#   link = "probit",
#   prior_intercept = normal(scale = 100),
#   prior_trt = normal(scale = 10),
#   prior_aux = flat(),
#   iter = 500,
#   chains = 2
# )
# 
# pasi_fit_fe_baseline <- nma(
#   pasi_net,
#   trt_effects = "fixed",
#   regression = ~ .mu,
#   link = "probit",
#   prior_intercept = normal(scale = 100),
#   prior_trt = normal(scale = 10),
#   prior_aux = flat(),
#   iter = niter
# )
# 
# results$pasi_fe <- nma_results(
#   pasi_fit_fe, 
#   base_dist = beta_dist_metaprop(pasi_ref, "random")
# )
# 
# results$pasi_fe_baseline <- nma_results(
#   pasi_fit_fe,
#   base_dist = beta_dist_metaprop(pasi_ref, "random"),
#   method = "baseline adjusted"
# )
# 
# 
# # Random effects
# pasi_fit_re <- nma(
#   pasi_net,
#   trt_effects = "random",
#   link = "probit",
#   prior_intercept = normal(scale = 100),
#   prior_trt = normal(scale = 10),
#   prior_aux = flat(),
#   iter = niter
# )
# 
# pasi_fit_re_baseline <- nma(
#   pasi_net,
#   trt_effects = "random",
#   regression = ~ .mu:.trt,
#   link = "probit",
#   prior_intercept = normal(scale = 100),
#   prior_trt = normal(scale = 10),
#   prior_aux = flat(),
#   iter = niter
# )
# 
# results$pasi_re <- nma_results(
#   pasi_fit_re, 
#   beta_dist_metaprop(pasi_ref, "random")
# )
# 
## DLQI response ---------------------------------------------------------------

dlqi_ref <- metaprop(
  event = dlqi_0_1,
  n = n,
  data = filter(data, drug == "Placebo", !is.na(dlqi_0_1)),
  sm = "PLOGIT",
  method = "Inverse",
  method.incr = "all",
  incr = 0.5
)

dlqi_net <- set_agd_arm(
  filter(data, !if_all(dlqi_0_1:dlqi_0, \(x) is.na(x))),
  study = ref_id,
  trt = drug,
  r =  multi(r0 = n,
             dlqi_0_1, dlqi_0,
             inclusive = TRUE,
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

results$dlqi_fe <- nma_results(
  dlqi_fit_fe, 
  beta_dist_metaprop(dlqi_ref, "random")
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

results$dlqi_re <- nma_results(
  dlqi_fit_re, 
  beta_dist_metaprop(dlqi_ref, "random")
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
  "pasi50", "pasi75", "pasi90", "pasi100", "sae", "disc_any", "disc_ae", 
  "serious_infection", "injection_site_rxn", "malignancy"
)

bin_fit_fe <- list()
bin_fit_re <- list()

for (i in 1:length(bin_outcomes)) {
  
  bin_net <- set_agd_arm(
    filter(data, !is.na(.data[[bin_outcomes[i]]])),
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
    # sm = "PLOGIT",
    # method = "Inverse",
    # method.incr = "all",
    # incr = 0.5,
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
      studylab = univar$ref_id,
      # sm = "PLOGIT",
      # method = "Inverse",
      # method.incr = "all",
      # incr = 0.5
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
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

exc <- c("Izokibep", "Mirikizumab", "Phototherapy")

results$pasi_fe <- NULL
results$pasi_re <- NULL
results$pasi_fe_baseline <- NULL

model_info <- bind_rows(results) |> 
  mutate(endpoint_group = if_else(likelihood == "multinomial" & 
                              str_detect(endpoint, "pasi|dlqi"),
                            substr(endpoint, 1, 4), endpoint)) |> 
  distinct(endpoint_group, type, likelihood, method, effects, dic) |> 
  arrange(type != "network", likelihood != "multinomial") |> 
  mutate(ma_id = row_number())

results_table <- bind_rows(results) |> 
  filter(comp_tx %notin% exc, ref_tx %notin% exc) |> 
  mutate(endpoint_group = if_else(likelihood == "multinomial" & 
                                    str_detect(endpoint, "pasi|dlqi"),
                                  substr(endpoint, 1, 4), endpoint)) |> 
  left_join(model_info, by = c("endpoint_group", "type", "likelihood", "method", "effects", "dic"),
            relationship = "many-to-one") |> 
  select(-c(type, likelihood, method, effects, dic, endpoint_group))

dbWriteTable(con, name = "ma_models", value = model_info, overwrite = TRUE)
dbWriteTable(con, name = "ma_results", value = results_table, overwrite = TRUE)

create_view_sql <- "
  CREATE VIEW v_meta_analysis AS
  SELECT *
  FROM ma_results
  LEFT JOIN ma_models ON ma_results.ma_id = ma_models.ma_id
"
dbExecute(con, "DROP VIEW IF EXISTS v_meta_analysis")
dbExecute(con, create_view_sql)

dbDisconnect(con)
save.image("R/meta-analyse/meta-analysis.RData")
