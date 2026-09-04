rm(list = ls())
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
theme_set(theme_classic())
j <- list()
j <- readRDS("R/meta-analyse/jags_fits.rds")

j$fe_fez_u_nc_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "fixed")
j$re_fez_u_nc_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "fixed")
j$fe_rezt_u_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "treatment")
j$re_rezt_u_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "treatment")
j$fe_rezi_u_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "trial")
j$re_rezi_u_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "trial")
j$fe_reza_u_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "arm")
j$re_reza_u_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "arm")

j$fe_fez_a_nc_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "fixed", baseline = "adjusted")
j$re_fez_a_nc_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "fixed", baseline = "adjusted")
j$fe_rezt_a_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "treatment", baseline = "adjusted")
j$re_rezt_a_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "treatment", baseline = "adjusted")
j$fe_rezi_a_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "trial", baseline = "adjusted")
j$re_rezi_a_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "trial", baseline = "adjusted")
j$fe_reza_a_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "arm", baseline = "adjusted")
j$re_reza_a_nc_con <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "arm", baseline = "adjusted")

j$fe_fez_u_c_con   <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "fixed", class = "exchangeable")
j$re_fez_u_c_con   <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "fixed", class = "exchangeable")
j$fe_rezt_u_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "treatment", class = "exchangeable")
j$re_rezt_u_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "treatment", class = "exchangeable")
j$fe_rezi_u_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "trial", class = "exchangeable")
j$re_rezi_u_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "trial", class = "exchangeable")
j$fe_reza_u_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "arm", class = "exchangeable")
j$re_reza_u_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "arm", class = "exchangeable")
j$fe_fez_a_c_con   <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "fixed", baseline = "adjusted", class = "exchangeable")
j$re_fez_a_c_con   <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "fixed", baseline = "adjusted", class = "exchangeable")
j$fe_rezt_a_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "treatment", baseline = "adjusted", class = "exchangeable")
j$re_rezt_a_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "treatment", baseline = "adjusted", class = "exchangeable")
j$fe_rezi_a_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "trial", baseline = "adjusted", class = "exchangeable")
j$re_rezi_a_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "trial", baseline = "adjusted", class = "exchangeable")
j$fe_reza_a_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "arm", baseline = "adjusted", class = "exchangeable")
j$re_reza_a_c_con  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "arm", baseline = "adjusted", class = "exchangeable")

j$fe_fez_u_nc_ume  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "fixed", consistency = "ume")
j$re_fez_u_nc_ume  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "fixed", consistency = "ume")
j$fe_rezt_u_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "treatment", consistency = "ume")
j$re_rezt_u_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "treatment", consistency = "ume")
j$fe_rezi_u_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "trial", consistency = "ume")
j$re_rezi_u_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "trial", consistency = "ume")
j$fe_reza_u_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "arm", consistency = "ume")
j$re_reza_u_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "arm", consistency = "ume")
j$fe_fez_a_nc_ume  <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "fixed", baseline = "adjusted", consistency = "ume")
j$re_fez_a_nc_ume  <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "fixed", baseline = "adjusted", consistency = "ume")
j$fe_rezt_a_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "treatment", baseline = "adjusted", consistency = "ume")
j$re_rezt_a_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "treatment", baseline = "adjusted", consistency = "ume")
j$fe_rezi_a_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "trial", baseline = "adjusted", consistency = "ume")
j$re_rezi_a_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "trial", baseline = "adjusted", consistency = "ume")
j$fe_reza_a_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "fixed", cutpoints = "arm", baseline = "adjusted", consistency = "ume")
j$re_reza_a_nc_ume <- pso_jags(pasi_jags, niter = 3000, effects = "random", cutpoints = "arm", baseline = "adjusted", consistency = "ume")

process_jags(j$re_rezi_u_nc_con)$dev_table

process_jags(j$re_rezi_u_nc_con)$dev_table |> 
  ggplot() +
  geom_density(aes(mean)) +
  geom_vline(xintercept = 1, linetype = 2)

devplot(j$re_rezi_u_nc_con, j$re_rezi_u_nc_ume)

devdev <- devplot(j$re_rezi_u_nc_con, j$re_rezi_u_nc_ume, "table")
devdev |> 
  filter(ref_id==112)


saveRDS(j, "R/meta-analyse/jags_fits.rds")


pasi |> 
  group_by(ref_id) |> 
  summarise(nt = n_distinct(drug)) |> 
  ungroup() |> 
  count(nt)


