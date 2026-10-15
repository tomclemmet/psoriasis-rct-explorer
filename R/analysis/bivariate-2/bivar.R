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
  arrange(ref_id, drug != "Placebo", drug) |> 
  mutate(.by = ref_id, arm_no = row_number(), .after = "ref_id")

# Extract DLQI data and assign treatment ids
dlqi <- dbReadTable(con, "v_dlqi") |> 
  filter(!is.na(dlqi_0_1)) |> 
  select(trial, ref_id, drug, n, dlqi_0_1, dlqi_0_1) |> 
  summarise(.by = c(ref_id, drug), n_dlqi = sum(n), dlqi_0_1 = sum(dlqi_0_1)) |> 
  arrange(ref_id, drug != "Placebo", drug) |> 
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

# Merge data and produce index lookup
bivar_data <- full_join(pasi, dlqi, by = c("ref_id", "arm_no", "drug")) |> 
  mutate(t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))))
t_lookup <- distinct(bivar_data, drug, t) |> arrange(t)
  
ggplot(bivar_data, aes(x = pasi75/n_pasi, y = dlqi_0_1/n_dlqi)) + 
  geom_point(aes(colour = drug, size = n_pasi)) +
  scale_colour_viridis_d()


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
  left_join(t_lookup, by = "drug") |> 
  mutate(
    .by = ref_id, 
    arm_no = arm_no,
    t = t,
    na2 = max(arm_no, na.rm = TRUE),
    o1 = if_else(!is.na(pasi75), 1, 2), # 1 for PASI 75, 2 for DLQI
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
  select(-c(ref_id, drug, t)) |> rename(arm = arm_no, study = s) |>
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
  nt = n_distinct(bivar_data$drug),
  na1 = row_lookup$na1,
  s = df2$s,
  t = select(df2, t1:t3) |> as.matrix(),
  na2 = df2$na2,
  o = select(df2, o1:o2) |> as.matrix(),
  oc = df3$oc,
  study = df3$study,
  arm = df3$arm,
  y = select(df3, y1:y2) |> as.matrix(),
  se = select(df3, se1:se2) |> as.matrix()
)
       
lapply(bivar_jags, head)

m <- jags(
  data = bivar_jags, parameters.to.save = c("d", "sd", "rhoB", "sig_t", "rho_t", "beta_t", "psi_t", "eta"),
  inits = NULL, model.file = "R/analysis/bivariate-2/bivar2d.jags", n.chains = 2, 
  n.iter = 20000, n.thin = 1
)
m

m_ind <- jags(
  data = bivar_jags, parameters.to.save = c("d", "sd"),
  inits = NULL, model.file = "R/analysis/bivariate-2/bivar2d-uni.jags", n.chains = 2, 
  n.iter = 20000, n.thin = 1
)

s <- m$BUGSoutput$summary |> 
  as_tibble(rownames = "param") |> 
  filter(str_starts(param, "d\\[")) |> 
  mutate(outcome = c("PASI 75", "DLQI")[as.numeric(str_sub(param, 3, 3))],
         t = as.numeric(str_extract(param, "(?<=,)[^]]+")),
         drug = if_else(outcome == "PASI 75", t_lookup$drug[t], t_lookup$drug[t]),
         .before = "param")

s_ind <- m_ind$BUGSoutput$summary |> 
  as_tibble(rownames = "param") |> 
  filter(str_starts(param, "d\\[")) |> 
  mutate(outcome = c("PASI 75", "DLQI")[as.numeric(str_sub(param, 3, 3))],
         t = as.numeric(str_extract(param, "(?<=,)[^]]+")),
         drug = if_else(outcome == "PASI 75", t_lookup$drug[t], t_lookup$drug[t]),
         .before = "param") |> 
  filter(drug %in% unique(dlqi$drug) | outcome == "PASI 75")

order <- s |> 
  filter(outcome == "PASI 75") |> 
  arrange(mean) |> 
  pull(drug) |> 
  unique()

bind_rows(s_ind, s, .id = "model") |> 
  mutate(drug = factor(drug, levels = order)) |> 
  ggplot(aes(x = mean, y = drug)) +
  geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = model), 
                  position = position_dodge(0.5, reverse = TRUE), 
                  size = 0.3, shape = 15) +
  facet_wrap(~ outcome, scale = "free_x") +
  scale_colour_discrete(labels = c("1" = "Univariate", "2" = "Bivariate")) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("output/bivar.png", height = 7, width = 12)

bind_rows(s_ind, s, .id = "model") |> 
  filter(drug != "Placebo") |> 
  arrange(outcome, drug, model) |> group_by(outcome, drug) |> 
  summarise(se_reduction = 1 - (last(sd)/first(sd)), .groups = "drop_last") |>
  summarise(mean_se_reduction = mean(se_reduction))

  
                  