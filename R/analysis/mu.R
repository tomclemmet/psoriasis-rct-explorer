source("R/utils/load-pasi.R")

rezi_u <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
)

rezi_a <- pso_jags(
  jags_inpits, niter = 200000, effects = random, cutpoints = "trial", 
  reg = "mu", reg_mean = -0.7, reg_scale = 1
)

rezt_u <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "treatment"
)

rezt_a <- pso_jags(
  jags_inputs, niter = 200000, effects = random, cutpoints = "treatment", 
  reg = "mu", reg_mean = -0.7, reg_scale = 1
)

out <- list(reziu = process_jags(rezi_u), rezia=process_jags(rezi_a), reztu = process_jags(rezt_u), rezt_a = process_jags(rezt_a))

saveRDS("output/mutest.rds")