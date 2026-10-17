rm(list = ls())
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R")

j <- list()
char <- "bio"
avg <- 0.25
sc <- 1
jags_inputs <- load_chars(char)

j$base <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 0.7), reg_scale = c(4, 1), 
  reg_effects = c("independent", "common"), reg_level = c("class", "common")
)
j$com <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", char), reg_mean = c(16, 0.7, avg), reg_scale = c(4, 1, sc), 
  reg_effects = c("independent", "common", "common"), 
  reg_level = c("class", "common", "common")
)
j$ex_class<- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", char), reg_mean = c(16, 0.7, avg), reg_scale = c(4, 1, sc), 
  reg_effects = c("independent", "common", "exchangeable"), 
  reg_level = c("class", "common", "class")
)
j$ind_class <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu", char), reg_mean = c(16, 0.7, avg), reg_scale = c(4, 1, sc), 
  reg_effects = c("independent", "common", "independent"), 
  reg_level = c("class", "common", "class")
)

saveRDS(lapply(j, process_jags), paste0("output/metareg/", char, "_mu.rds"))