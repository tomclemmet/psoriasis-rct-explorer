rm(list = ls())
library(meta)
library(ggplot2)
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
theme_set(theme_classic())
bayesplot::color_scheme_set("viridis")
j <- list()
ume <- list()
# j <- readRDS("R/meta-analyse/jags_fits.rds")
niter <- 10000

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

ume <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", consistency = "ume", baseline = "adjusted")
j$re_rezi_a_nc
ume

devplot(j$re_rezi_a_nc, ume, "Consistency", "Inconsistency")

# STEP 3: Check class effects model

class <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", consistency = "ume", baseline = "adjusted", class = "exchangeable")
j$re_rezi_a_nc
class

# STEP 4: Sensitivity analysis
s <- list()

source("R/meta-analyse/wide_format.R")
s$base <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")

bios <- gen_pasi_jags(c("Placebo", "Adalimumab", "Bimekizumab", "Brodalumab", "Certolizumab",
                        "Etanercept", "Guselkumab", "Infliximab", "Ixekizumab", "Risankizumab", "Secukinumab", "Sonelokimab",
                        "Tildrakizumab", "Ustekinumab", "Xeligekimab", "Netakimab", "Mirikizumab"))
s$bio <- pso_jags(bios, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")

source("R/meta-analyse/wide_format.R")
approved <- gen_pasi_jags(setdiff(pasi_drugs, c(
  "Icotrokinra", "Mirikizumab", "Netakimab", "Orismilast", "Roflumilast", 
  "Phototherapy", "Sonelokimab", "Tofacitinib", "Xeligekimab")
  ))
s$app <- pso_jags(approved, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")

source("R/meta-analyse/wide_format.R")
s$low_rob <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")

forests(s$base, s$bio, s$app, s$low_rob) +
  scale_colour_viridis_d(labels = c("1" = "base case", "2" = "biologics", "3" = "approved", "4" = "low RoB"), end = 0.8)
ggsave("output/sa_forest.png", height = 18, width = 16, units = "cm")

# Compare predictions of cutpoints models
j$re_fez_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed", baseline = "adjusted")
j$re_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted")
j$re_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted")
j$re_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm", baseline = "adjusted")

forests(j)
lapply(j, \(x) x)

j[[4]]$trace |> mcmc_trace(pars = "B")

umet <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted", consistency = "ume")
umei <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial", baseline = "adjusted", consistency = "ume")


devplot(j$re_rezt_a_nc, umet, "Consistency", "Inconsistency")

devplot(j$re_rezi_a_nc, umei, "Consistency", "Inconsistency")
