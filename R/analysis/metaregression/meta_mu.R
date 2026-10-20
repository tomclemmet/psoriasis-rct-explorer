rm(list = ls())
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R")

j <- list()
char <- "age"
avg <- 45
sc <- 10
jags_inputs <- load_chars(char)


j$age <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "age"), reg_mean = c(16, 0.7, 45), reg_scale = c(4, 1, 10), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$sex <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "sex"), reg_mean = c(16, 0.7, 0.7), reg_scale = c(4, 1, 1), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$ethnicity <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "white"), reg_mean = c(16, 0.7, 0.9), reg_scale = c(4, 1, 1), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$weight <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "weight"), reg_mean = c(16, 0.7, 90), reg_scale = c(4, 1, 10), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$psa <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "psa"), reg_mean = c(16, 0.7, 0.25), reg_scale = c(4, 1, 1), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$duration <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "duration"), reg_mean = c(16, 0.7, 15), reg_scale = c(4, 1, 5), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$pasi <- pso_jags(
  jags_inputs, niter = 60000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", "pasi"), reg_mean = c(16, 0.7, 15), reg_scale = c(4, 1, 5), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)


saveRDS(lapply(j, process_jags), paste0("output/metareg/meta_mu.rds"))