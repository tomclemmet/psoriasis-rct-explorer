rm(list = ls())
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
theme_set(theme_classic())
j <- list()
ume <- list()
# j <- readRDS("R/meta-analyse/jags_fits.rds")
niter <- 4000

# j$fe_fez_u_nc  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed")
# j$re_fez_u_nc  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed")
# j$fe_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment")
# j$re_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment")
# j$fe_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial")
j$re_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial")
# j$fe_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm")
# j$re_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm")

# j$fe_fez_a_nc  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed", baseline = "adjusted")
# j$re_fez_a_nc  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", baseline = "adjusted")
# j$fe_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment", baseline = "adjusted")
# j$re_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted")
# j$fe_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial", baseline = "adjusted")
# j$re_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")
# j$fe_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm", baseline = "adjusted")
# j$re_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", baseline = "adjusted")
# 
# j$fe_fez_u_c   <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed", class = "exchangeable")
# j$re_fez_u_c   <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", class = "exchangeable")
# j$fe_rezt_u_c  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment", class = "exchangeable")
# j$re_rezt_u_c  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", class = "exchangeable")
# j$fe_rezi_u_c  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial", class = "exchangeable")
# j$re_rezi_u_c  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", class = "exchangeable")
# j$fe_reza_u_c  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm", class = "exchangeable")
# j$re_reza_u_c  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", class = "exchangeable")
# j$fe_fez_a_c   <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed", baseline = "adjusted", class = "exchangeable")
# j$re_fez_a_c   <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", baseline = "adjusted", class = "exchangeable")
# j$fe_rezt_a_c  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment", baseline = "adjusted", class = "exchangeable")
# j$re_rezt_a_c  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted", class = "exchangeable")
# j$fe_rezi_a_c  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial", baseline = "adjusted", class = "exchangeable")
# j$re_rezi_a_c  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted", class = "exchangeable")
# j$fe_reza_a_c  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm", baseline = "adjusted", class = "exchangeable")
# j$re_reza_a_c  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", baseline = "adjusted", class = "exchangeable")

# ume$fe_fez_u_nc  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed", consistency = "ume")
# ume$re_fez_u_nc  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", consistency = "ume")
# ume$fe_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment", consistency = "ume")
# ume$re_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", consistency = "ume")
# ume$fe_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial", consistency = "ume")
ume$re_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", consistency = "ume")
# ume$fe_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm", consistency = "ume")
# ume$re_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", consistency = "ume")
# ume$fe_fez_a_nc  <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "fixed", baseline = "adjusted", consistency = "ume")
# ume$re_fez_a_nc  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", baseline = "adjusted", consistency = "ume")
# ume$fe_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "treatment", baseline = "adjusted", consistency = "ume")
# ume$re_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted", consistency = "ume")
# ume$fe_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "trial", baseline = "adjusted", consistency = "ume")
# ume$re_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted", consistency = "ume")
# ume$fe_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed", cutpoints = "arm", baseline = "adjusted", consistency = "ume")
# ume$re_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", baseline = "adjusted", consistency = "ume")

# summ <- compare_jags(j)
# summ
# rez_rezi_u_nc is the best in terms of dic and 2nd best in terms of resdev

process_jags(j$re_rezi_u_nc)$dev_table |> filter(ref_id==85)
process_jags(ume$re_rezi_u_nc)
# UME model has slightly lower DIC, though this is driven by pV to a large extent
# Worth investigating signs of inconsistency

devplot(j$re_rezi_u_nc, ume$re_rezi_u_nc)
# Several points are significantly away from the line

devplot(j$re_rezi_u_nc, ume$re_rezi_u_nc, "table") |> 
  filter(diff > 0.5) |> mutate(.by = ref_id, tdiff = mean(diff)) |> arrange(desc(tdiff))
# 13 / 958 data points have over 0.5 points of inconsistency
# Fixed error around CCA/NRI in Allah-farani which may help reich
# Removed khalid for an unapproved dose
# Need to investigate BROvsGUS loops more, note this includes Reich
# Pruned dose violations globally which may help


saveRDS(j, "R/meta-analyse/jags_fits.rds")




