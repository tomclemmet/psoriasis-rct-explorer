library(mvnma)

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")


pasi <- dbReadTable(con, "v_pasi") |> 
  filter(!is.na(pasi75)) |> 
  summarise(.by = c(trial, ref_id, drug), n = sum(n), pasi75 = sum(pasi75))

pw1 <- pairwise(treat = drug, event = pasi75, n = n, studlab = ref_id, data = pasi, sm = "RD", reference.group = "Placebo")

dlqi <- dbReadTable(con, "v_dlqi") |> 
  filter(!is.na(dlqi_0_1)) |> 
  summarise(.by = c(trial, ref_id, drug), n = sum(n), dlqi_0_1 = sum(dlqi_0_1))

pw2 <- pairwise(treat = drug, event = dlqi_0_1, n = n, studlab = ref_id, data = dlqi, sm = "OR", reference.group = "Placebo")

m <- mvnma(pw1, pw2, reference.group = "Placebo", outclab = c("PASI 75", "DLQI 0/1"),
           n.iter = 5000, n.chains = 2, lower.rho = 0.4, upper.rho = 0.6)

m0 <- netmeta(pw2)

m
as.mcmc(m) |> bayesplot::mcmc_trace(pars = "rho")
m0
