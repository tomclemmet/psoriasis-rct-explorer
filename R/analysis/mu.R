source("R/utils/load-pasi.R")
source("R/utils/ord-jags-nma.R")

rezi_u <- pso_jags(
  jags_inputs, niter = 5000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4,
  reg_effects = c("independent"), reg_level = c("class")
)

rezi_a <- pso_jags(
  jags_inputs, niter = 200000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "mu"), reg_mean = c(16, -0.7), reg_scale = c(4, 1),
  reg_effects = c("independent", "common"), reg_level = c("class", "common")
)

rezt_u <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "treatment",
  reg = "timepoint", reg_mean = 16, reg_scale = 4,
  reg_effects = c("independent"), reg_level = c("class")
)

rezt_a <- pso_jags(
  jags_inputs, niter = 200000, effects = "random", cutpoints = "treatment", 
  reg = c("timepoint", "mu"), reg_mean = c(16, -0.7), reg_scale = c(4, 1),
  reg_effects = c("independent", "common"), reg_level = c("class", "common")
)

out <- list(rezi_u = process_jags(rezi_u), rezi_a=process_jags(rezi_a), 
            reztu = process_jags(rezt_u), rezt_a = process_jags(rezt_a))

saveRDS(out, "output/mu_test.rds")
