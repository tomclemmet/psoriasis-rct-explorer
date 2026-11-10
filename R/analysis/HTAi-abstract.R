rm(list = ls())
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-pasi.R")

settings <- expand_grid(
  effects = c("fixed", "random"),
  cutpoints = c("fixed", "treatment", "trial",  "arm"),
  class = c("independent", "exchangeable"),
  reg = list("timepoint", c("timepoint", "mu")),
  reg_mean = list(c(16), c(16, -0.7)),
  reg_scale = list(c(4), c(4, 1)),
  reg_effects = list(c("independent"), c("independent", "common")),
  reg_level = list(c("class"), c("class", "treatment"))
) |> 
  filter(if_all(starts_with("reg_"), \(x) lengths(x) == lengths(reg)))

names <- expand_grid(
  effects = c("fe", "re"),
  cutpoints = c("fez", "rezt", "rezi",  "rezi"),
  class = c("ind", "cl"),
  reg = c("unadj", "mu"),
) |> 
  mutate(name = paste(effects, cutpoints, class, reg, sep = "_")) |> 
  pull(name)

m <- list()
names(m) <- names
for (i in 1:nrow(settings)) {
  mu <- settings[[i, "reg"]][[1]][2] == "mu"
  m[[i]] <- pso_jags(
    jags_inputs,
    niter = 100000,
    effects = settings[[i, "effects"]],
    cutpoints = settings[[i, "cutpoints"]],
    class = settings[[i, "class"]],
    reg = settings[[i, "reg"]][[1]],
    reg_mean = settings[[i, "reg_mean"]][[1]],
    reg_scale = settings[[i, "reg_scale"]][[1]],
    reg_effects = settings[[i, "reg_effects"]][[1]],
    reg_level = settings[[i, "reg_level"]][[1]]
  )
}

p <- lapply(m, jags_process)
saveRDS(p, "output/htai.rds")

p <- readRDS("C:/Users/p12916tc/Documents/Data/NMA/saved-results/htai.rds")

compare_jags(p)









