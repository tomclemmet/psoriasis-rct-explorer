rm(list = ls())
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-pasi.R")

settings <- expand_grid(
  effects = "fixed",
  cutpoints = c("fixed", "treatment", "trial",  "arm"),
  class = c("independent", "exchangeable"),
  reg = list("timepoint", c("timepoint", "mu")),
  reg_mean = list(c(16), c(16, -0.7)),
  reg_scale = list(c(4), c(4, 1)),
  reg_effects = list(c("independent"), c("independent", "common")),
  reg_level = list(c("class"), c("class", "treatment"))
) |> bind_rows(
    expand_grid(
      effects = "random",
      cutpoints = c("fixed", "treatment", "trial",  "arm"),
      class = c("independent", "exchangeable"),
      reg = list("timepoint", c("timepoint", "mu")),
      reg_mean = list(c(16), c(16, -1)),
      reg_scale = list(c(4), c(4, 1)),
      reg_effects = list(c("independent"), c("independent", "common")),
      reg_level = list(c("class"), c("class", "treatment"))
    )
  ) |> 
  filter(if_all(starts_with("reg_"), \(x) lengths(x) == lengths(reg)))

names <- expand_grid(
  effects = c("fe", "re"),
  cutpoints = c("fez", "rezt", "rezi",  "reza"),
  class = c("ind", "cl"),
  reg = c("unadj", "mu"),
) |> 
  mutate(name = paste(effects, cutpoints, class, reg, sep = "_")) |> 
  pull(name)

m <- list()
for (i in 1:nrow(settings)) {
  mu <- settings[[i, "reg"]][[1]][2] == "mu"
  m[[i]] <- pso_jags(
    jags_inputs,
    # filename = paste0("JAGS/", names[i], ".jags"),
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
names(p) <- names
saveRDS(p, "output/htai-pasi.rds")

# p <- readRDS("C:/Users/p12916tc/Documents/Data/NMA/saved-results/htai.rds")
# 
# p2 <- lapply(p, \(x) {
#   x$summary <- x$summary |> 
#     mutate(
#       label = case_when(
#         param == "beta_timepoint[1]" ~ "beta_timepoint[placebo]",
#         param == "beta_timepoint[2]" ~ "beta_timepoint[bio_il1223]",
#         param == "beta_timepoint[3]" ~ "beta_timepoint[bio_il17]",
#         param == "beta_timepoint[4]" ~ "beta_timepoint[bio_il23]",
#         param == "beta_timepoint[5]" ~ "beta_timepoint[bio_tnf]",
#         param == "beta_timepoint[6]" ~ "beta_timepoint[non-biologic]",
#         .default = label
#       )
#     ) |> 
#     filter(!str_starts(param, "beta_mu"), param != "B_timepoint")
#   
#   # x$results <- x$results |> 
#   #   mutate(
#   #     label = case_when(
#   #       param == "beta_timepoint[1]" ~ "beta_timepoint[placebo]",
#   #       param == "beta_timepoint[2]" ~ "beta_timepoint[bio_il1223]",
#   #       param == "beta_timepoint[3]" ~ "beta_timepoint[bio_il17]",
#   #       param == "beta_timepoint[4]" ~ "beta_timepoint[bio_il23]",
#   #       param == "beta_timepoint[5]" ~ "beta_timepoint[bio_tnf]",
#   #       param == "beta_timepoint[6]" ~ "beta_timepoint[non-biologic]",
#   #       .default = "label"
#   #     )
#   #   ) |> 
#   #   filter(!str_starts(param, "beta_mu"), param != "B_timepoint")
#   x
# })
# compare_jags(p2) |> write.csv("output/htai.csv")
# 
# 
# forests(p2)
# ggsave("output/htai.png", height = 30, width = 10)
 
source("R/utils/load-dlqi.R")

settings <- expand_grid(
  effects = "fixed",
  cutpoints = c("fixed", "treatment", "trial",  "arm"),
  class = c("independent", "exchangeable"),
  reg = list("timepoint", c("timepoint", "mu")),
  reg_mean = list(c(16), c(16, -1.3)),
  reg_scale = list(c(4), c(4, 1)),
  reg_effects = list(c("independent"), c("independent", "common")),
  reg_level = list(c("class"), c("class", "treatment"))
) |> bind_rows(
  expand_grid(
    effects = "random",
    cutpoints = c("fixed", "treatment", "trial",  "arm"),
    class = c("independent", "exchangeable"),
    reg = list("timepoint", c("timepoint", "mu")),
    reg_mean = list(c(16), c(16, -1.3)),
    reg_scale = list(c(4), c(4, 1)),
    reg_effects = list(c("independent"), c("independent", "common")),
    reg_level = list(c("class"), c("class", "treatment"))
  )
) |> 
  filter(if_all(starts_with("reg_"), \(x) lengths(x) == lengths(reg)))

d <- list()
for (i in 1:nrow(settings)) {
  mu <- settings[[i, "reg"]][[1]][2] == "mu"
  d[[i]] <- pso_jags(
    jags_inputs,
    # filename = paste0("JAGS/", names[i], ".jags"),
    niter = 100,
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

dp <- lapply(d, jags_process)
names(dp) <- names
saveRDS(dp, "output/htai-dlqi.rds")

