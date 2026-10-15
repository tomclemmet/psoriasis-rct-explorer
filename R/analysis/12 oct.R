source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R") 
source("R/utils/load-pasi.R")


j <- list()
jp <- list()

j$unadjusted <- pso_jags(
  pasi_jags, niter = 15000, effects = "random", cutpoints = "trial"
)

j$common <- pso_jags(
  pasi_jags, niter = 15000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 1, reg_effects = "common"
)

j$ex <- pso_jags(
  pasi_jags, niter = 15000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 1, 
  reg_effects = "exchangeable", reg_level = "class"
)

j$ind <- pso_jags(
  pasi_jags, niter = 15000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 1, 
  reg_effects = "independent", reg_level = "class"
)

j$unadjusted_mu <- pso_jags(
  pasi_jags, niter = 50000, effects = "random", cutpoints = "trial",
  reg = "mu", reg_mean = 0.7, reg_scale = 1
)

j$common_mu <- pso_jags(
  pasi_jags, niter = 50000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 0.7), reg_scale = c(1, 1), 
  reg_effects = c("common", "common")
)

j$ext_mu <- pso_jags(
  pasi_jags, niter = 50000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 0.7), reg_scale = c(1, 1), 
  reg_effects = c("exchangeable", "common"), reg_level = c("class", "common")
)

j$ind_mu <- pso_jags(
  pasi_jags, niter = 50000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 0.7), reg_scale = c(1, 1), 
  reg_effects = c("independent", "common"), reg_level = c("class", "common")
)
