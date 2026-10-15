# Bivariate NMA code based on approach in TSD20 appendix D.2 (p.116)
options(digits = 3)
rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(R2jags)
library(ggplot2)
theme_set(theme_classic())

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

# Extract PASI data and assign treatment ids
pasi <- dbReadTable(con, "v_pasi") |> 
  filter(!is.na(pasi75)) |> 
  select(trial, ref_id, drug, n, pasi75) |> 
  summarise(.by = c(ref_id, drug), n_pasi = sum(n), pasi75 = sum(pasi75)) |> 
  mutate(t_pasi = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))), .after = "drug") |> 
  arrange(ref_id, t_pasi) |> 
  mutate(.by = ref_id, arm_no = row_number(), .after = "ref_id")

# Extract DLQI data and assign treatment ids
dlqi <- dbReadTable(con, "v_dlqi") |> 
  filter(!is.na(dlqi_0_1)) |> 
  select(trial, ref_id, drug, n, dlqi_0_1, dlqi_0_1) |> 
  summarise(.by = c(ref_id, drug), n_dlqi = sum(n), dlqi_0_1 = sum(dlqi_0_1)) |> 
  mutate(t_dlqi = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))), .after = "drug") |> 
  arrange(ref_id, t_dlqi) |> 
  mutate(.by = ref_id, arm_no = row_number(), .after = "ref_id")
dbDisconnect(con)

# dlqi_jags <- list(
#   ns = n_distinct(dlqi$ref_id),
#   nt = max(dlqi$t_dlqi),
#   na = summarise(dlqi, .by = ref_id, n = n())$n,
#   t = select(dlqi, ref_id, arm_no, t_dlqi) |> 
#     pivot_wider(names_from = arm_no, values_from = t_dlqi, names_prefix = "a") |> 
#     select(-ref_id) |> as.matrix(),
#   r = select(dlqi, ref_id, arm_no, dlqi_0_1) |> 
#     pivot_wider(names_from = arm_no, values_from = dlqi_0_1, names_prefix = "a") |> 
#     select(-ref_id) |> as.matrix(),
#   n = select(dlqi, ref_id, arm_no, n_dlqi) |> 
#     pivot_wider(names_from = arm_no, values_from = n_dlqi, names_prefix = "a") |> 
#     select(-ref_id) |> as.matrix()
# )

# Produce index lookups for both outcomes
p_index <- distinct(pasi, drug, t_pasi) |> arrange(t_pasi)
d_index <- distinct(dlqi, drug, t_dlqi) |> arrange(t_dlqi)

# Merge PASI and DLQI data
bivar_data <- full_join(pasi, dlqi, by = c("ref_id", "arm_no", "drug"))

bivar_data |> 
  mutate(drug = toupper(str_sub(drug, 1, 3))) |> 
  ggplot(aes(x = pasi75/n_pasi, y = dlqi_0_1/n_dlqi)) + 
  geom_point(aes(colour = drug, size = n_pasi)) +
  scale_colour_viridis_d() +
  guides(size = "none")

bivar_data |> 
  mutate(y1 = pasi75/n_pasi, y2 = dlqi_0_1/n_dlqi) |> 
  filter(!is.na(y1), !is.na(y2)) |> 
  summarise(c = cor(y1, y2))

row_lookup <- bivar_data |> 
  arrange(
    !is.na(pasi75) & !is.na(dlqi_0_1), 
    is.na(pasi75),
    ref_id
  ) |> 
  summarise(.by = ref_id, na1 = n()) |> 
  mutate(s = row_number())


# Produce datafile 2, containing treatment indices and outcome indicators
df2 <- bind_rows(pasi, dlqi) |> 
  mutate(
    .by = ref_id, 
    arm_no = arm_no,
    t = coalesce(t_pasi, t_dlqi), 
    na2 = max(arm_no, na.rm = TRUE),
    o1 = if_else(!is.na(t_pasi), 1, 2), # 1 for PASi 75, 2 for DLQI
    o2 = o1, # Usually equals o1
    .keep = "none"
  ) |> 
  pivot_wider(names_from = arm_no, values_from = t, names_prefix = "t") |>
  bind_rows(data.frame( # Add in rows for trials with PASI but no DLQI
    ref_id = setdiff(pasi$ref_id, dlqi$ref_id),
    o1 = 2,
    o2 = 0
  )) |> 
  bind_rows(data.frame( # Add in rows for trials with DLQI but no PASI
    ref_id = setdiff(dlqi$ref_id, pasi$ref_id),
    o1 = 1,
    o2 = 0
  )) |>
  left_join(row_lookup, by = "ref_id") |> 
  select(-ref_id) |> 
  arrange(s, o1) |> fill(c(na2, t1:t3)) |> relocate(s, t1:t3)


df3 <- bivar_data |> 
  left_join(row_lookup, by = "ref_id") |> 
  select(-c(ref_id, drug, t_pasi, t_dlqi)) |> rename(arm = arm_no, study = s) |>
  mutate(across(n_pasi:dlqi_0_1, \(x) {x + 0.5}),
         y1 = log(pasi75 / (n_pasi - pasi75)), y2 = log(dlqi_0_1 / (n_dlqi - dlqi_0_1)),
         se1 = sqrt(1/pasi75 + 1/(n_pasi - pasi75)), se2 = sqrt(1/dlqi_0_1 + 1/(n_dlqi - dlqi_0_1)),
         v = NA, oc = if_else(is.na(y1), 2, 1)) |>
  arrange(study, arm)

bivar_jags <- list(
  ns = n_distinct(bivar_data$ref_id),
  no = 2,
  nobs1 = sum(is.na(bivar_data$pasi75) | is.na(bivar_data$dlqi_0_1)),
  nobs2 = sum(!is.na(bivar_data$pasi75) & !is.na(bivar_data$dlqi_0_1)),
  nt.total = c(n_distinct(pasi$drug), n_distinct(dlqi$drug)),
  na1 = row_lookup$na1,
  s = df2$s,
  t = select(df2, t1:t3) |> as.matrix(),
  na2 = df2$na2,
  o = select(df2, o1:o2) |> as.matrix(),
  oc = df3$oc,
  study = df3$study,
  arm = df3$arm,
  y = select(df3, y1:y2) |> as.matrix(),
  se = select(df3, se1:se2) |> as.matrix()#,
  # v = df3$v
)
       
lapply(bivar_jags, head)

m <- jags(
  data = bivar_jags, parameters.to.save = c("d", "sd", "rhoB", "rhoW"),
  inits = NULL, model.file = "R/analysis/bivariate/bivar.jags", n.chains = 2, 
  n.iter = 10000, n.thin = 1
)
m

m_ind <- jags(
  data = bivar_jags, parameters.to.save = c("d", "sd"),
  inits = NULL, model.file = "R/analysis/bivariate/bivar-ind.jags", n.chains = 2, 
  n.iter = 5000, n.thin = 1
)

s <- m$BUGSoutput$summary |> 
  as_tibble(rownames = "param") |> 
  filter(str_starts(param, "d\\[")) |> 
  mutate(outcome = c("PASI 75", "DLQI")[as.numeric(str_sub(param, 3, 3))],
         t = as.numeric(str_extract(param, "(?<=,)[^]]+")),
         drug = if_else(outcome == "PASI 75", p_index$drug[t], d_index$drug[t]),
         .before = "param")

s_ind <- m_ind$BUGSoutput$summary |> 
  as_tibble(rownames = "param") |> 
  filter(str_starts(param, "d\\[")) |> 
  mutate(outcome = c("PASI 75", "DLQI")[as.numeric(str_sub(param, 3, 3))],
         t = as.numeric(str_extract(param, "(?<=,)[^]]+")),
         drug = if_else(outcome == "PASI 75", p_index$drug[t], d_index$drug[t]),
         .before = "param")

order <- s_ind |> 
  filter(outcome == "PASI 75") |> 
  arrange(mean) |> 
  pull(drug)

bind_rows(s, s_ind, .id = "model") |> 
  group_by(drug) |> 
  filter(!is.na(drug), any(outcome == "DLQI")) |> 
  mutate(drug = factor(drug, levels = order)) |> 
  ggplot(aes(x = mean, y = drug)) +
  geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = model), 
                  position = position_dodge(0.5), show.legend = F) +
  facet_wrap(~ outcome, scale = "free_x") +
  scale_colour_discrete(labels = c("1" = "Bivariate", "2" = "Univariate")) +
  theme_bw()

bivar_data |> 
  group_by(ref_id) |> filter(any(drug == "Brodalumab"), any(drug == "Placebo"), 
                             drug %in% c("Brodalumab", "Placebo"))
                  