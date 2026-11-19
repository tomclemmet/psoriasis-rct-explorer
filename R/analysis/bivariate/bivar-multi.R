# Multivariate bivariate NMA code
options(digits = 3)
rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(R2jags)
library(ggplot2)
theme_set(theme_classic())
classes <- read.csv("R/trt_class.csv")

con <- dbConnect(RSQLite::SQLite(), "pso-trial-explorer/psoriasis-rcts.sqlite")

source("R/utils/load-pasi.R")
pasi <- jags_inputs
pasi_id <- id_lookup
pasi_ref <- tibble(p = 1:pasi$ns, id = pasi_wide$ref_id)
names(pasi) <- paste0("pasi_", names(pasi))
source("R/utils/load-dlqi.R")
dlqi <- jags_inputs
dlqi_id <- id_lookup
dlqi_ref <- tibble(d = 1:dlqi$ns, id = dlqi_wide$ref_id)
names(dlqi) <- paste0("dlqi_", names(dlqi))

comb <- c(pasi, dlqi)
comb$tt <- pasi_id |> 
  select(drug, t) |>
  right_join(select(dlqi_id, drug, t), by = "drug") |> 
  pull(t.x)
# NOTE: Only works if all DLQI trials report PASI
# This object can be used to map between DLQI and PASI trials. For items 
# 1:dlqi_ns the number reports to corresponding the PASI trial row number,
# the remaining items are the remaining PASI trial row numbers where only 
# PASI is reported
comb$ref_map <- left_join(pasi_ref, dlqi_ref, by = "id") |> 
  relocate(id) |> arrange(is.na(d)) |> 
  pull(p)


m <- list()

# m$fe_ind <- jags(
#   data = comb, parameters.to.save = c("pasi_d", "dlqi_d"),
#   inits = NULL, model.file = "R/analysis/bivariate-2/bivar-multi-fe-ind.jags", n.chains = 2, 
#   n.iter = 5000, n.thin = 1
# )
# m$fe_cov <- jags(
#   data = comb, parameters.to.save = c("pasi_d", "dlqi_d", "rho_t"),
#   inits = NULL, model.file = "R/analysis/bivariate-2/bivar-multi-fe-cov.jags", n.chains = 2, 
#   n.iter = 15000, n.thin = 1
# )
m$re_ind <- jags(
  data = comb, parameters.to.save = c("pasi_d", "dlqi_d", "pasi_sd", "dlqi_sd"),
  inits = NULL, model.file = "JAGS/bivariate/bivar-multi-re-ind.jags", n.chains = 2, 
  n.iter = 5000, n.thin = 1
)
m$re_cov <- jags(
  data = comb, parameters.to.save = c("pasi_d", "dlqi_d", "pasi_sd", "dlqi_sd", "rhoB", "rho_t"),
  inits = NULL, model.file = "JAGS/bivariate/bivar-multi-re-cov.jags", n.chains = 2, 
  n.iter = 15000, n.thin = 1
)
# m$fe_ind
# m$fe_cov
m$re_ind
m$re_cov



summ <- function(mod) {
  mod$BUGSoutput$summary |> 
    as_tibble(rownames = "param") |> 
    filter(str_detect(param, "d\\[")) |> 
    mutate(outcome = str_sub(param, 1, 4),
           t = as.numeric(str_sub(param, 8, -2)),
           drug = if_else(outcome == "PASI 75", pasi_id$drug[t], pasi_id$drug[t]),
           .before = "param")
}

s <- lapply(m, summ)


order <- s$re_ind |> 
  filter(outcome=="pasi") |> 
  arrange(mean) |> 
  pull(drug) |> 
  unique()

bind_rows(s, .id = "model") |> 
  mutate(drug = factor(drug, levels = order),
         model = factor(model, levels = c("re_ind", "re_cov"))) |> 
  # mutate(across(c(mean, `2.5%`, `97.5%`), exp)) |> 
  ggplot(aes(x = mean, y = drug)) +
  geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = model), 
                  position = position_dodge(0.5, reverse = TRUE), 
                  size = 0.25, shape = 15) +
  facet_wrap(~ outcome, scale = "free_x") +
  theme_bw() +
  theme(legend.position = "right") +
  labs(title = "Bivariate and univariate NMA results", x = "Probit response vs placebo", y = "")
ggsave("output/bivar_multi.png", height = 15, width = 30, units = "cm")











