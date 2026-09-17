rm(list = ls())
library(meta)
library(ggplot2)
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
theme_set(theme_classic())
j <- list()
ume <- list()
# j <- readRDS("R/meta-analyse/jags_fits.rds")
niter <- 4000

# STEP 1: Compare random and fixed effects models, with and without baseline adjustment
j$fe_fez_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed")
j$re_fez_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed")
j$fe_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment")
j$re_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment")
j$fe_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial")
j$re_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial")
j$fe_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm")
j$re_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm")
j$fe_fez_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed", baseline = "adjusted")
j$re_fez_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", baseline = "adjusted")
j$fe_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment", baseline = "adjusted")
j$re_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted")
j$fe_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial", baseline = "adjusted")
j$re_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")
j$fe_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm", baseline = "adjusted")
j$re_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", baseline = "adjusted")
saveRDS(j, "R/meta-analyse/jags_fits.rds")

compare_jags(j) |> write.csv("results1.csv")

# STEP 2: Check consistency using UME model
ume$re_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", consistency = "ume", baseline = "adjusted")

m <- j$re_rezi_a_nc
u <- ume$re_rezi_a_nc
m
u
forest(m)

m$dev_table |> arrange(desc(mean))

devplot(m, u, "Consistency", "Inconsistency")

devplot(m, u, output = "table") |> 
  filter(diff > 0.5) |> group_by(ref_id) |> mutate(mdiff = max(diff)) |> arrange(mdiff) |> View()

distinct(pasi_condensed, ref_id, trial) |> write.csv("temp.csv")

