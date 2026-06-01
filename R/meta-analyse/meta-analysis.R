rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)
library(meta)

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")
dbListTables(con)

query <- "
SELECT 
    p.ref_id, p.arm_no, p.drug, p.timepoint, p.timepoint_unit, p.n,
    p.pasi50, p.pasi75, p.pasi90, p.pasi100, p.abs_pasi_change_mean, 
    p.abs_pasi_change_sd,
    d.dlqi_0_1, d.dlqi_0, d.abs_dlqi_change_mean, d.abs_dlqi_change_sd,
    s.sae, s.disc_any, s.disc_ae, s.serious_infection, s.injection_site_rxn,
    s.malignancy
FROM v_pasi p
LEFT JOIN v_dlqi d      ON p.ref_id = d.ref_id   AND p.arm_no = d.arm_no
LEFT JOIN v_safety s    ON p.ref_id = s.ref_id   AND p.arm_no = s.arm_no
"

data <- dbGetQuery(con, query) |> filter(!is.na(drug))
results <- data.frame(NULL)

# PASI Response NMA ===========================================================
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
  iter = 100
)

pasi_ref_fe <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(pasi_fit_fe, pars = "mu"))
)
pasi_rd_fe <- predict(
    pasi_fit_fe, type = "response",
    baseline = distr(qnorm, pasi_ref_fe$TE.fixed, pasi_ref_fe$seTE.fixed)
  )$sims |> 
  posterior::as_draws_df() |>
  mutate(across(contains("pasi50"), \(x) x - `pred[Placebo, pasi50]`)) |> 
  mutate(across(contains("pasi75"), \(x) x - `pred[Placebo, pasi75]`)) |> 
  mutate(across(contains("pasi90"), \(x) x - `pred[Placebo, pasi90]`)) |> 
  mutate(across(contains("pasi100"), \(x) x - `pred[Placebo, pasi100]`)) |>
  pivot_longer(everything(), names_to = "param", values_to = "trace") |>
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])")) |> 
  filter(drug != "Placebo")
pasi_rate_fe <- predict(
  pasi_fit_fe, type = "response",
  baseline = distr(qnorm, pasi_ref_fe$TE.fixed, pasi_ref_fe$seTE.fixed)
)$sims |> 
  posterior::as_draws_df() |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |>  
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])"))

results <- results |> bind_rows(data.frame(
    endpoint = pasi_rd_fe$endpoint,
    type = "network",
    effects = "fixed",
    ref_tx = "Placebo",
    comp_tx = pasi_rd_fe$drug,
    measure = "rd",
    mean = pasi_rd_fe$mean,
    lower = pasi_rd_fe$lower,
    upper = pasi_rd_fe$upper
  ), data.frame(
    endpoint = pasi_rate_fe$endpoint,
    type = "network",
    effects = "fixed",
    ref_tx = NA,
    comp_tx = pasi_rate_fe$drug,
    measure = "rate",
    mean = pasi_rate_fe$mean,
    lower = pasi_rate_fe$lower,
    upper = pasi_rate_fe$upper
  ))


pasi_fit_re <- nma(
  pasi_net,
  trt_effects = "random",
  link = "probit",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  prior_aux = flat(),
  iter = 100
)
pasi_ref_re <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(pasi_fit_fe, pars = "mu"))
)
pasi_rd_re <- predict(
  pasi_fit_re, type = "response",
  baseline = distr(qnorm, pasi_ref_re$TE.random, pasi_ref_re$seTE.random)
)$sims |> 
  posterior::as_draws_df() |>
  mutate(across(contains("pasi50"), \(x) x - `pred[Placebo, pasi50]`)) |> 
  mutate(across(contains("pasi75"), \(x) x - `pred[Placebo, pasi75]`)) |> 
  mutate(across(contains("pasi90"), \(x) x - `pred[Placebo, pasi90]`)) |> 
  mutate(across(contains("pasi100"), \(x) x - `pred[Placebo, pasi100]`)) |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])")) |> 
  filter(drug != "Placebo")
pasi_rate_re <- predict(
  pasi_fit_re, type = "response",
  baseline = distr(qnorm, pasi_ref_re$TE.random, pasi_ref_re$seTE.random)
)$sims |> 
  posterior::as_draws_df() |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])"))

results <- results |> bind_rows(data.frame(
  endpoint = pasi_rd_re$endpoint,
  type = "network",
  effects = "random",
  ref_tx = "Placebo",
  comp_tx = pasi_rd_re$drug,
  measure = "rd",
  mean = pasi_rd_re$mean,
  lower = pasi_rd_re$lower,
  upper = pasi_rd_re$upper
), data.frame(
  endpoint = pasi_rate_re$endpoint,
  type = "network",
  effects = "random",
  ref_tx = NA,
  comp_tx = pasi_rate_re$drug,
  measure = "rate",
  mean = pasi_rate_re$mean,
  lower = pasi_rate_re$lower,
  upper = pasi_rate_re$upper
))

# DLQI response NMA ===========================================================
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
  iter = 100
)
dlqi_ref_fe <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(dlqi_fit_fe, pars = "mu"))
)
dlqi_rd_fe <- predict(
  dlqi_fit_fe, type = "response",
  baseline = distr(qnorm, dlqi_ref_fe$TE.fixed, dlqi_ref_fe$seTE.fixed)
)$sims |> 
  posterior::as_draws_df() |> 
  mutate(across(contains("dlqi50"), \(x) x - `pred[Placebo, dlqi50]`)) |> 
  mutate(across(contains("dlqi75"), \(x) x - `pred[Placebo, dlqi75]`)) |> 
  mutate(across(contains("dlqi90"), \(x) x - `pred[Placebo, dlqi90]`)) |> 
  mutate(across(contains("dlqi100"), \(x) x - `pred[Placebo, dlqi100]`)) |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])")) |> 
  filter(drug != "Placebo")
dlqi_rate_fe <- predict(
  dlqi_fit_fe, type = "response",
  baseline = distr(qnorm, dlqi_ref_fe$TE.fixed, dlqi_ref_fe$seTE.fixed)
)$sims |> 
  posterior::as_draws_df() |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])"))

results <- results |> bind_rows(data.frame(
  endpoint = dlqi_rd_fe$endpoint,
  type = "network",
  effects = "fixed",
  ref_tx = "Placebo",
  comp_tx = dlqi_rd_fe$drug,
  measure = "rd",
  mean = dlqi_rd_fe$mean,
  lower = dlqi_rd_fe$lower,
  upper = dlqi_rd_fe$upper
), data.frame(
  endpoint = dlqi_rate_fe$endpoint,
  type = "network",
  effects = "fixed",
  ref_tx = NA,
  comp_tx = dlqi_rate_fe$drug,
  measure = "rate",
  mean = dlqi_rate_fe$mean,
  lower = dlqi_rate_fe$lower,
  upper = dlqi_rate_fe$upper
))

dlqi_fit_re <- nma(
  dlqi_net,
  trt_effects = "random",
  link = "probit",
  prior_intercept = normal(scale = 100),
  prior_trt = normal(scale = 10),
  prior_aux = flat(),
  iter = 100
)
dlqi_ref_re <- metagen(
  TE = mean, 
  seTE = sd, 
  data = as.data.frame(summary(dlqi_fit_fe, pars = "mu"))
)
dlqi_rd_re <- predict(
  dlqi_fit_re, type = "response",
  baseline = distr(qnorm, dlqi_ref_re$TE.random, dlqi_ref_re$seTE.random)
)$sims |> 
  posterior::as_draws_df() |> 
  mutate(across(contains("dlqi50"), \(x) x - `pred[Placebo, dlqi50]`)) |> 
  mutate(across(contains("dlqi75"), \(x) x - `pred[Placebo, dlqi75]`)) |> 
  mutate(across(contains("dlqi90"), \(x) x - `pred[Placebo, dlqi90]`)) |> 
  mutate(across(contains("dlqi100"), \(x) x - `pred[Placebo, dlqi100]`)) |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])")) |> 
  filter(drug != "Placebo")
dlqi_rate_re <- predict(
  dlqi_fit_re, type = "response",
  baseline = distr(qnorm, dlqi_ref_re$TE.random, dlqi_ref_re$seTE.random)
)$sims |> 
  posterior::as_draws_df() |> 
  pivot_longer(everything(), names_to = "param", values_to = "trace") |> 
  summarise(.by = param, mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)) |> 
  filter(substr(param, 1, 1) != ".") |> 
  mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
         endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])"))

results <- results |> bind_rows(data.frame(
  endpoint = dlqi_rd_re$endpoint,
  type = "network",
  effects = "random",
  ref_tx = "Placebo",
  comp_tx = dlqi_rd_re$drug,
  measure = "rd",
  mean = dlqi_rd_re$mean,
  lower = dlqi_rd_re$lower,
  upper = dlqi_rd_re$upper
), data.frame(
  endpoint = dlqi_rate_re$endpoint,
  type = "network",
  effects = "random",
  ref_tx = NA,
  comp_tx = dlqi_rate_re$drug,
  measure = "rate",
  mean = dlqi_rate_re$mean,
  lower = dlqi_rate_re$lower,
  upper = dlqi_rate_re$upper
))

# Pairwise Meta-Analyses ======================================================
outcomes <- cbind(
  c("pasi50", "pasi75", "pasi90", "pasi100", "abs_pasi_change_mean",
    "dlqi_0_1", "dlqi_0", "sae", "disc_any", "disc_ae", "serious_infection", 
    "injection_site_rxn", "malignancy"),
  c("binary", "binary", "binary", "binary", "continuous", "binary", "binary", 
    "binary", "binary", "binary", "binary", "binary", "binary")
)

comparisons <- data |>
  distinct(ref_id, drug) |> arrange(ref_id, drug) |> 
  inner_join(select(data, ref_id, drug), by = "ref_id",
             relationship = "many-to-many") |> 
  filter(drug.x < drug.y) |> 
  distinct(drug.x, drug.y)

for (i in 1:nrow(outcomes)) {
  for (j in 1:nrow(comparisons)) {
    comp_data <- data |> 
      mutate(
        k = .data[[outcomes[i, 1]]], 
        lab = if_else(drug == comparisons[[j, 1]], 1, 2)
      ) |> 
      select(ref_id, lab, n, k) |> 
      filter(!is.na(k)) |> 
      group_by(ref_id) |> 
      filter(
        !is.na(lab),
        any(lab == 1) & any(lab == 2),
      ) |> ungroup() |> 
      summarise(.by = c(ref_id, lab), n = sum(n), k = sum(k)) |> 
      tidyr::pivot_wider(names_from = lab, values_from = c(n, k))
    
    if(nrow(comp_data) == 0 || n_distinct(comp_data$ref_id) <= 1) next
      
    if (outcomes[i, 2] == "binary") {
      fit <- metabin(
        event.e = comp_data$k_1, n.e = comp_data$n_1,event.c = comp_data$k_2, 
        n.c = comp_data$n_2, sm = "RD"
      )
      abs_fit_1 <- metaprop(event = comp_data$k_1, n = comp_data$n_1)
      abs_fit_2 <- metaprop(event = comp_data$k_2, n = comp_data$n_2)
    }
    results <- results |> bind_rows(data.frame(
      endpoint = outcomes[i, 1],
      type = "pairwise",
      effects = "fixed",
      ref_tx = comparisons[[j, 1]],
      comp_tx = comparisons[[j, 2]],
      outcome = "rd",
      mean = fit$TE.common,
      lower = fit$lower.common,
      upper = fit$upper.common
    ), data.frame(
      endpoint = outcomes[i, 1],
      type = "pairwise",
      effects = "random",
      ref_tx = comparisons[[j, 1]],
      comp_tx = comparisons[[j, 2]],
      outcome = "rd",
      mean = fit$TE.random,
      lower = fit$lower.random,
      upper = fit$upper.random
    ), data.frame(
      endpoint = outcomes[i, 1],
      type = "univariate",
      effects = "fixed",
      ref_tx = NA,
      comp_tx = comparisons[[j, 1]],
      outcome = "rate",
      mean = plogis(abs_fit_1$TE.fixed),
      lower = plogis(abs_fit_1$lower.fixed),
      upper = plogis(abs_fit_1$upper.fixed)
    ), data.frame(
      endpoint = outcomes[i, 1],
      type = "univariate",
      effects = "random",
      ref_tx = NA,
      comp_tx = comparisons[[j, 1]],
      outcome = "rate",
      mean = plogis(abs_fit_1$TE.random),
      lower = plogis(abs_fit_1$lower.random),
      upper = plogis(abs_fit_1$upper.random)
    ), data.frame(
      endpoint = outcomes[i, 1],
      type = "univariate",
      effects = "fixed",
      ref_tx = NA,
      comp_tx = comparisons[[j, 2]],
      outcome = "rate",
      mean = plogis(abs_fit_2$TE.fixed),
      lower = plogis(abs_fit_2$lower.fixed),
      upper = plogis(abs_fit_2$upper.fixed)
    ), data.frame(
      endpoint = outcomes[i, 1],
      type = "univariate",
      effects = "random",
      ref_tx = NA,
      comp_tx = comparisons[[j, 2]],
      outcome = "rate",
      mean = plogis(abs_fit_2$TE.random),
      lower = plogis(abs_fit_2$lower.random),
      upper = plogis(abs_fit_2$upper.random)
    ))
    message(".", appendLF = FALSE)
  }
}


# if (spec$kind == "binary")
#   meta::metabin(event.e = pt$event_a, n.e = pt$n_a,
#                 event.c = pt$event_b, n.c = pt$n_b,
#                 studlab = pt$trial, sm = "RR", method = "Inverse",
#                 method.tau = "REML", common = TRUE, random = TRUE,
#                 warn = FALSE)
# else
#   meta::metacont(n.e = pt$n_a, mean.e = pt$mean_a, sd.e = pt$sd_a,
#                  n.c = pt$n_b, mean.c = pt$mean_b, sd.c = pt$sd_b,
#                  studlab = pt$trial, sm = "MD", method.tau = "REML",
#                  common = TRUE, random = TRUE, warn = FALSE)


dbWriteTable(con, name = "meta_analysis", value = results, overwrite = TRUE)

dbDisconnect()