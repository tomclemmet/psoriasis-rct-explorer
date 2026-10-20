rm(list = ls())
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R")

j <- list()
jc <- list()
n_iter = 50000

jags_inputs <- load_chars("age")
j$age <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$age <- process_jags(j$age)
jags_inputs <- load_chars("male")
j$male <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$male <- process_jags(j$male)
jags_inputs <- load_chars("white")
j$white <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$white <- process_jags(j$white)
jags_inputs <- load_chars("weight")
j$weight <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$weight <- process_jags(j$weight)
jags_inputs <- load_chars("psa")
j$psa <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$psa <- process_jags(j$psa)
jags_inputs <- load_chars("duration")
j$duration <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$duration <- process_jags(j$duration)
jags_inputs <- load_chars("pasi")
j$pasi <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$pasi <- process_jags(j$pasi)
jags_inputs <- load_chars("bio")
j$bio <- pso_jags(
  jags_inputs, niter = n_iter, effects = "random", cutpoints = "trial",
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("independent"), 
  reg_level = c("class")
)
jc$bio <- process_jags(j$bio)

saveRDS(jc, "output/metareg/unadjusted.rds")