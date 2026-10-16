source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-pasi.R")

j <- list()

j$base <- pso_jags(
  jags_inputs, niter = 10000, effects = "random", cutpoints = "trial"
)
j$com <- pso_jags(
  jags_inputs, niter = 10000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4
)
j$ex_class<- pso_jags(
  jags_inputs, niter = 10000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "exchangeable",
  reg_level = "class"
)
j$ind_class <- pso_jags(
  jags_inputs, niter = 10000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "independent",
  reg_level = "class"
)


check_fits(j)
write.csv(compare_jags(j), "output/metareg/timepoint.csv")
forests(j)
ggsave("output/metareg/timepoint.png", height = 10, width = 7)

j_mu <- list()

j_mu$base <- pso_jags(
  jags_inputs, niter = 30000, effects = "random", cutpoints = "trial",
  reg = "mu", reg_mean = 0.7, reg_scale = 1
)
j_mu$com <- pso_jags(
  jags_inputs, niter = 30000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 1), reg_scale = c(4, 1)
)
j_mu$ex_class<- pso_jags(
  jags_inputs, niter = 30000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 1), reg_scale = c(4, 1), 
  reg_effects = c("exchangeable", "common"), reg_level = c("class", "common")
)
j_mu$ind_class <- pso_jags(
  jags_inputs, niter = 30000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "mu"), reg_mean = c(16, 1), reg_scale = c(4, 1), 
  reg_effects = c("independent", "common"), reg_level = c("class", "common")
)


check_fits(j_mu)
write.csv(compare_jags(j_mu), "output/metareg/timepoint_mu.csv")
forests(j_mu)
ggsave("output/metareg/timepoint_mu.png", height = 10, width = 7)
