rm(list = ls())
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
theme_set(theme_classic())
j <- list()
ume <- list()
# j <- readRDS("R/meta-analyse/jags_fits.rds")
niter <- 4000

distinct(pasi_condensed, ref_id, trial, pop_res) |> 
  filter(!is.na(pop_res) & pop_res != "null") |> 
  arrange(pop_res)

j$re_rezi_u_nc <- pso_jags(gen_pasi_jags(), niter = niter, effects = "random", cutpoints = "trial")
ume$re_rezi_u_nc <- pso_jags(gen_pasi_jags(), niter = niter, effects = "random", cutpoints = "trial", consistency = "ume")

traceplot(j$re_rezi_u_nc, varname = "beta")

m <- process_jags(j$re_rezi_u_nc)
u <- process_jags(ume$re_rezi_u_nc)
m
u
forest(u)

devplot(m, u, xlab="Consistency", ylab="Inconsistency")

devplot(m, u, "table") |> 
  filter(diff > 0.5) |> group_by(ref_id) |> mutate(mdiff = max(diff)) |> arrange(desc(diff)) |> View()#|> distinct(trial) |> left_join(pasi_condensed)


t1 <- "Risankizumab"
t2 <- "Secukinumab"

direct_data <- gen_pasi_jags(c(t1, t2))
if (nrow(pasi_wide) == 1) {
  direct_data <- list(
    nc = pasi_wide$nc,
    C = direct_data$C |> as.vector(),
    r = rbind(
      select(pasi_wide, starts_with("a1r")) |> as.vector(),
      select(pasi_wide, starts_with("a2r")) |> as.vector()
    ),
    n = rbind(
      select(pasi_wide, starts_with("a1n")) |> as.vector(),
      select(pasi_wide, starts_with("a2n")) |> as.vector()
    )
  )
  direct <- jags(
    data = direct_data,
    parameters.to.save = "d", inits = NULL, model.file = "JAGS/single_study.jags", 
    n.chains = 2, n.iter = 100000, n.burnin = 50000
  )
  dir <- posterior::as_draws_df(direct$BUGSoutput$sims.array)$d
} else {
  direct <- pso_jags(direct_data, effects = "random", cutpoints = "trial")
  dir <- posterior::as_draws_df(direct$BUGSoutput$sims.array)$`d[2]`
}

indirect_data <- gen_pasi_jags(c(t1, t2), direct = FALSE)
indirect <- pso_jags(indirect_data, effects = "random", cutpoints = "trial")

p1 <- paste0("d[", match(t1, pasi_drugs), "]")
p2 <- paste0("d[", match(t2, pasi_drugs), "]")

indirect_trace <- posterior::as_draws_df(indirect$BUGSoutput$sims.array)

data.frame(t1 = indirect_trace[[p1]], t2 = indirect_trace[[p2]]) |> 
  mutate(indir = t2 - t1, dir = dir) |> 
  select(-t1, -t2) |> 
  pivot_longer(everything(), names_to = "lab", values_to = "d") |> 
  ggplot() +
  geom_density(aes(d, colour = lab, fill = lab), alpha = 0.5)

gen_pasi_jags()

# The indirect evidence is pulling RIS up and SEC down
# Align on timepoints to fix immerge
# Check ECLIPSE for SEc > GUS
# Check reason to exclude Immprint

pasi_condensed |> 
  select(-pasi_high_rob, -pop_res) |> 
  group_by(ref_id) |> 
  filter(any(drug == "Secukinumab")) |> 
  ungroup() |> 
  filter(drug != "Secukinumab") |> 
  mutate(control = "Secukinumab") |> 
  left_join(select(pasi_condensed, ref_id, drug, n, pasi50, pasi75, pasi90, pasi100), by = c("ref_id", "control" = "drug")) |> 
  group_by(drug) |> 
  arrange(drug) |> 
  summarise(
    n = sum(n.x + n.y),
    p50 = if_else(n() == 1, 1/sum((pasi50.x/n.x) / (pasi50.y/n.y)), 1/exp(metabin(pasi50.x, n.x, pasi50.y, n.y, data = pick(everything()))$TE.random)),
    p75 = if_else(n() == 1, 1/sum((pasi75.x/n.x) / (pasi75.y/n.y)), 1/exp(metabin(pasi75.x, n.x, pasi75.y, n.y, data = pick(everything()))$TE.random)),
    p90 = if_else(n() == 1, 1/sum((pasi90.x/n.x) / (pasi90.y/n.y)), 1/exp(metabin(pasi90.x, n.x, pasi90.y, n.y, data = pick(everything()))$TE.random)),
    p100 = if_else(n() == 1, 1/sum((pasi100.x/n.x) / (pasi100.y/n.y)), 1/exp(metabin(pasi100.x, n.x, pasi100.y, n.y, data = pick(everything()))$TE.random))
  ) |> 
  arrange(desc(n))

pasi_condensed |> 
  select(-pasi_high_rob, -pop_res) |> 
  group_by(ref_id) |> 
  filter(any(drug == "Risankizumab")) |> 
  ungroup() |> 
  filter(drug != "Risankizumab") |> 
  mutate(control = "Risankizumab") |> 
  left_join(select(pasi_condensed, ref_id, drug, n, pasi50, pasi75, pasi90, pasi100), by = c("ref_id", "control" = "drug")) |> 
  group_by(drug) |> 
  arrange(drug) |> 
  summarise(
    n = sum(n.x + n.y),
    p50 = if_else(n() == 1, 1/sum((pasi50.x/n.x) / (pasi50.y/n.y)), 1/exp(metabin(pasi50.x, n.x, pasi50.y, n.y, data = pick(everything()))$TE.random)),
    p75 = if_else(n() == 1, 1/sum((pasi75.x/n.x) / (pasi75.y/n.y)), 1/exp(metabin(pasi75.x, n.x, pasi75.y, n.y, data = pick(everything()))$TE.random)),
    p90 = if_else(n() == 1, 1/sum((pasi90.x/n.x) / (pasi90.y/n.y)), 1/exp(metabin(pasi90.x, n.x, pasi90.y, n.y, data = pick(everything()))$TE.random)),
    p100 = if_else(n() == 1, 1/sum((pasi100.x/n.x) / (pasi100.y/n.y)), 1/exp(metabin(pasi100.x, n.x, pasi100.y, n.y, data = pick(everything()))$TE.random))
  ) |> 
  arrange(desc(n))

pasi_condensed |> 
  select(-pasi_high_rob, -pop_res) |> 
  group_by(ref_id, trial) |> 
  filter(any(drug == "Secukinumab") & any(drug == "Placebo"), drug %in% c("Secukinumab", "Placebo")) |> 
  arrange(drug, .by_group = TRUE) |> 
  summarise(rd50 = last(pasi50)/last(n) - first(pasi50)/first(n),
            rd75 = last(pasi75)/last(n) - first(pasi75)/first(n),
            rd90 = last(pasi90)/last(n) - first(pasi90)/first(n),
            rd100 = last(pasi100)/last(n) - first(pasi100)/first(n))

sec <- pasi_condensed |> 
  select(-pasi_high_rob, -pop_res) |> 
  group_by(ref_id, trial) |> 
  filter(any(drug == "Secukinumab") & any(drug == "Placebo"), drug %in% c("Secukinumab", "Placebo")) |> 
  arrange(drug, .by_group = TRUE) |> 
  mutate(group = last(drug),
         p050 = pasi50/n,
         p075 = pasi75/n,
         p090 = pasi90/n,
         p100 = pasi100/n) |> 
  select(-(pasi50:pasi100)) |> 
  pivot_longer(p050:p100, names_to = "response", values_to = "p")

ris <- pasi_condensed |> 
  select(-pasi_high_rob, -pop_res) |> 
  group_by(ref_id, trial) |> 
  filter(any(drug == "Risankizumab") & any(drug == "Placebo"), drug %in% c("Risankizumab", "Placebo")) |> 
  arrange(drug, .by_group = TRUE) |> 
  mutate(group = last(drug),
         p050 = pasi50/n,
         p075 = pasi75/n,
         p090 = pasi90/n,
         p100 = pasi100/n) |> 
  select(-(pasi50:pasi100)) |> 
  pivot_longer(p050:p100, names_to = "response", values_to = "p")

ggplot(bind_rows(sec, ris)) +
  geom_point(aes(y = p, x = group, colour = drug), position = position_jitter(width = 0.2, height = 0, seed = 123)) +
  facet_wrap(~ response, scale = "free_x") +
  scale_colour_viridis_d()
