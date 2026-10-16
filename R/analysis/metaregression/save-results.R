source("R/utils/plots.R")

m <- readRDS("R/saved-results/timepoint_mu.rds")


check_fits(m)

# m$ind_class$summary <- m$ind_class$summary |> filter(param!="B_timepoint")
# m$com$summary <- m$com$summary |> filter(!str_starts(param, "beta_mu"))
# m$base$summary <- m$base$summary |> filter(!str_starts(param, "beta_mu"))
# m$ex_class$summary <- m$ex_class$summary |> filter(!str_starts(param, "beta_mu"))
# m$ind_class$summary <- m$ind_class$summary |> filter(!str_starts(param, "beta_mu"))


write.csv(compare_jags(m), "output/metareg/timepoint_mu.csv")
forests(m)
ggsave("output/metareg/timepoint_mu.png", height = 10, width = 7)
