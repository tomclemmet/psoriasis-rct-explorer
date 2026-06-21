rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)

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

bin_outcomes <- c(
  "pasi50", "pasi75", "pasi90", "pasi100",
  "dlqi_0_1", "dlqi_0", 
  "sae", "disc_any", "disc_ae", "serious_infection", "injection_site_rxn", 
  "malignancy"
)

# Binary outcomes ==============================================================

for (i in 1:length(bin_outcomes)) {
  for (j in 1:nrow(comparisons)) {
    tx <- comparisons[[j, 1]]
    ref <- comparisons[[j, 2]]
    outcome <- bin_outcomes[i]
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

    if(nrow(comp_data) == 0) next
    
    results[[paste(outcome, tx, ref)]] <- comp_data |> mutate(
      endpoint = outcome,
      comp_tx = tx,
      ref_tx = ref,
      p_tx = k_tx / n_tx,
      p_ref = k_ref / n_ref,
      p_tx_se = sqrt(p_tx * (1 - p_tx) / n_tx),
      p_ref_se = sqrt(p_ref * (1 - p_ref) / n_ref),
      measure = "rd",
      mean = p_tx - p_ref,
      lower = pmax(-1, mean - 1.96 * sqrt(p_tx_se^2 + p_ref_se^2)),
      upper = pmin(1, mean + 1.96 * sqrt(p_tx_se^2 + p_ref_se^2))
    ) |> 
      select(-starts_with("p_"))
  }
}

for (i in 1:length(bin_outcomes)) {
  for (k in 1:length(drugs)) {
    results[[paste(bin_outcomes[i], drugs[k])]] <- data |> 
      filter(
        !is.na(.data[[bin_outcomes[i]]]),
        drug == drugs[k]
      ) |> 
      select(ref_id, drug, n, contains(bin_outcomes[i])) |> 
      summarise(.by = c(ref_id, drug), n = sum(n), k = sum(.data[[bin_outcomes[i]]])) |> 
      mutate(
        mean = k / n,
        lower = pmax(0, mean - 1.96 * sqrt(mean * (1 - mean) / n)),
        upper = pmin(1, mean + 1.96 * sqrt(mean * (1 - mean) / n))
      ) |> 
      rename(n_tx = n, k_tx = k, comp_tx = drug) |> 
      mutate(endpoint = bin_outcomes[i], measure = "rate")
      
  }
}

# Continuous outcomes ==========================================================

## Absolute PASI ---------------------------------------------------------------

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
  filter(!is.na(abs_pasi_change_mean) & !is.na(abs_pasi_change_sd)) |> 
  rename(comp_tx = drug) |> 
  summarise(.by = c(ref_id, comp_tx), 
            mean = mean(abs_pasi_change_mean), sd = mean(abs_pasi_change_sd), # Not quite accurate
            n = sum(n)) |> 
  mutate(
    measure = "cfb",
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) |> 
  select(ref_id, n, mean, sd, lower, upper, comp_tx, measure)

results$`abs_pasi_change cfb` <- abs_pasi_data |> 
  mutate(endpoint = "abs_pasi_change") |> 
  rename(n_tx = n, sd_tx = sd)

for (j in 1:nrow(comparisons)) {
  tx <- comparisons[[j, 1]]
  ref <- comparisons[[j, 2]]
  comp_data <- abs_pasi_data |> 
    group_by(ref_id) |> 
    filter(any(comp_tx == tx) & any(comp_tx == ref),
           comp_tx %in% c(tx, ref)) |> 
    ungroup() |> 
    mutate(label = if_else(comp_tx == tx, "tx", "ref")) |> 
    select(-c(lower, upper, measure)) |> 
    pivot_wider(names_from = label, values_from = c(comp_tx, n, mean, sd))
  
  if(nrow(comp_data) == 0) next
  
  results[[paste("abs_pasi_change", tx, ref)]] <- comp_data |> 
    mutate(measure = "diff_cfb", endpoint = "abs_pasi_change",
           mean = mean_tx - mean_ref, se = sqrt((sd_ref^2 + sd_tx^2) / (n_ref + n_tx)),
           lower = mean - 1.96 * se, upper = mean + 1.96 * se) |> 
    rename(comp_tx = comp_tx_tx, ref_tx = comp_tx_ref) |> select(-se)
}

## Absolute DLQI ---------------------------------------------------------------

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
  filter(!is.na(abs_dlqi_change_mean) & !is.na(abs_dlqi_change_sd)) |> 
  rename(comp_tx = drug) |> 
  summarise(.by = c(ref_id, comp_tx), 
            mean = mean(abs_dlqi_change_mean), sd = mean(abs_dlqi_change_sd), # Not quite accurate
            n = sum(n)) |> 
  mutate(
    measure = "cfb",
    lower = mean - 1.96 * sd / sqrt(n),
    upper = mean + 1.96 * sd / sqrt(n)
  ) |> 
  select(ref_id, n, mean, sd, lower, upper, comp_tx, measure)

results$`abs_dlqi_change cfb` <- abs_dlqi_data |> 
  mutate(endpoint = "abs_dlqi_change") |> 
  rename(n_tx = n, sd_tx = sd)

for (j in 1:nrow(comparisons)) {
  tx <- comparisons[[j, 1]]
  ref <- comparisons[[j, 2]]
  comp_data <- abs_dlqi_data |> 
    group_by(ref_id) |> 
    filter(any(comp_tx == tx) & any(comp_tx == ref),
           comp_tx %in% c(tx, ref)) |> 
    ungroup() |> 
    mutate(label = if_else(comp_tx == tx, "tx", "ref")) |> 
    select(-c(lower, upper, measure)) |> 
    pivot_wider(names_from = label, values_from = c(comp_tx, n, mean, sd))
  
  if(nrow(comp_data) == 0) next
  
  results[[paste("abs_dlqi_change", tx, ref)]] <- comp_data |> 
    mutate(measure = "diff_cfb", endpoint = "abs_pasi_change",
           mean = mean_tx - mean_ref, se = sqrt((sd_ref^2 + sd_tx^2) / (n_ref + n_tx)),
           lower = mean - 1.96 * se, upper = mean + 1.96 * se) |> 
    rename(comp_tx = comp_tx_tx, ref_tx = comp_tx_ref) |> select(-se)
}

results_table <- bind_rows(results)

dbWriteTable(con, name = "trial_estimates", value = results_table, overwrite = TRUE)

create_view_sql <- "
  CREATE VIEW v_trial_estimates AS
  SELECT *
  FROM trial_estimates
"
dbExecute(con, create_view_sql)

dbDisconnect(con)
