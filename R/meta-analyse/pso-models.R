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

m <- process_jags(j$re_rezi_u_nc)
u <- process_jags(ume$re_rezi_u_nc)
m
u

m$results |> arrange(desc(Rhat)) |> head(20)
u$results |> arrange(desc(Rhat)) |> head(20)

m$dev_table |> arrange(desc(mean)) |> head(20)
u$dev_table |> arrange(desc(mean)) |> head(20)

devplot(m, u, xlab="Consistency", ylab="Inconsistency")

devplot(m, u, "table") |> 
  filter(diff > 0.5) |> group_by(ref_id) |> mutate(mdiff = max(diff)) |> arrange(desc(diff)) #|> distinct(trial) |> left_join(pasi_condensed)


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

pasi_id_lookup

posterior::as_draws_df(indirect$BUGSoutput$sims.array) |> 
  mutate(indir = `d[21]` - `d[19]`, dir = dir, .keep = "none") |> 
  pivot_longer(everything(), names_to = "lab", values_to = "d") |> 
  ggplot() +
  geom_density(aes(d, colour = lab, fill = lab), alpha = 0.5)

# The indirect evidence is pulling RIS up and SEC down

