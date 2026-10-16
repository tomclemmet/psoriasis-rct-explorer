rm(list = ls())
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R")

j <- list()
char <- "pasi"
avg <- 15
sc <- 4
jags_inputs <- load_chars(char)

j$base <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "independent",
  reg_level = "class"
)
j$com <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", char), reg_mean = c(16, avg), reg_scale = c(4, sc), 
  reg_effects = c("independent", "common"), reg_level = c("class", "common")
)
j$ex_class<- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, avg), reg_scale = c(4, sc), 
  reg_effects = c("independent", "exchangeable"), reg_level = c("class", "class")
)
j$ind_class <- pso_jags(
  jags_inputs, niter = 20000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", char), reg_mean = c(16, avg), reg_scale = c(4, sc), 
  reg_effects = c("independent", "independent"), reg_level = c("class", "class")
)

saveRDS(lapply(j, process_jags), paste0("output/metareg/", char, ".rds"))