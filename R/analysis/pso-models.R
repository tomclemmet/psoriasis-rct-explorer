rm(list = ls())
library(ggplot2)
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
theme_set(theme_classic())
bayesplot::color_scheme_set("viridis")
j <- list()
ume <- list()
class <- list()
j <- readRDS("R/meta-analyse/jags_fits.rds")
ume <- readRDS("R/meta-analyse/jags_fits.rds")
class <- readRDS("R/meta-analyse/jags_fits.rds")
niter <- 10000


# STEP 1: Compare model fits ----------------------------------------------

j$fe_fez_u_nc  <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "fixed")
j$re_fez_u_nc  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed")
j$fe_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "treatment")
j$re_rezt_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment")
j$fe_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "trial")
j$re_rezi_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial")
j$fe_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "arm")
j$re_reza_u_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm")
j$fe_fez_a_nc  <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "fixed",     baseline = "adjusted")
j$re_fez_a_nc  <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "fixed",     baseline = "adjusted")
j$fe_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "treatment", baseline = "adjusted")
j$re_rezt_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "treatment", baseline = "adjusted")
j$fe_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "trial",     baseline = "adjusted")
j$re_rezi_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "trial",     baseline = "adjusted")
j$fe_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "fixed",  cutpoints = "arm",       baseline = "adjusted")
j$re_reza_a_nc <- pso_jags(pasi_jags, niter = niter, effects = "random", cutpoints = "arm",       baseline = "adjusted")
saveRDS(j, "R/meta-analyse/jags_fits.rds")
lapply(j, \(x) x$DIC)
compare_jags(j) |> write.csv("output/results1.csv")


forests(j$fe_fez_a_nc, j$re_rezt_a_nc, j$re_rezi_a_nc)
forests(j$fe_fez_a_nc, j$re_rezt_a_nc, j$re_rezi_a_nc, prob = TRUE)

# STEP 2: Check consistency -----------------------------------------------

# Treatment-level cutpoints
ume$t <- pso_jags(pasi_jags, niter = 10000, effects = "random", cutpoints = "treatment", baseline = "adjusted", consistency = "ume")
ume$t
j$re_rezt_a_nc
devplot(j$re_rezt_a_nc, ume$t, "Consistency", "Inconsistency")
ggsave("output/devdev_t.png", height = 5, width = 5)
devplot(j$re_rezt_a_nc, ume$t, output="table") |> 
  filter(diff > 0.5) |> arrange(-diff)

ume$i <- pso_jags(pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", baseline = "adjusted", consistency = "ume")
ume$i
j$re_rezi_a_nc
devplot(j$re_rezi_a_nc, ume$i, "Consistency", "Inconsistency")
ggsave("output/devdev_i.png", height = 5, width = 5)
devplot(j$re_rezi_a_nc, ume$i, output="table") |> 
  filter(diff > 0.5) |> arrange(-diff)

saveRDS(ume, "R/meta-analyse/ume_fits.rds")

# STEP 3: Check class effects ---------------------------------------------

class$rezt <- pso_jags(pasi_jags, niter = 10000, effects = "random", cutpoints = "treatment", baseline = "adjusted", class = "exchangeable")
class$rezt
j$re_rezt_a_nc
forests(class$rezt, j$re_rezt_a_nc)
forests(class$rezt, j$re_rezt_a_nc, prob = TRUE)

class$rezi <- pso_jags(pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", baseline = "adjusted", class = "exchangeable")
class$rezi
j$re_rezi_a_nc
forests(class$rezi, j$re_rezi_a_nc)
forests(class$rezi, j$re_rezi_a_nc, prob = TRUE)

class$fe <- pso_jags(pasi_jags, niter = 10000, effects = "fixed", cutpoints = "fixed", baseline = "adjusted", class = "exchangeable")
class$fe
j$fe_fez_a_nc
forests(class$fe, j$fe_fez_a_nc)
forests(class$fe, j$fe_fez_a_nc, prob = TRUE)

# Class effects don't meaningfully change estimates and DIC increases
saveRDS(ume, "R/meta-analyse/class_fits.rds")

# STEP 4: Sensitivity analysis --------------------------------------------

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

