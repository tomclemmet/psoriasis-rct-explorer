source("R/utils/plots.R")
source("R/utils/load-chars.R")

chars <- c("age", "male", "white", "weight", "psa", "duration", "pasi", "bio")

unadjusted <- readRDS("R/saved-results/metareg/unadjusted.rds")
mu <- readRDS("R/saved-results/metareg/mu.rds")
metamu <- readRDS("R/saved-results/metareg/meta_mu.rds")

for (i in 1:length(unadjusted)) {
  mu[[i]]$summary <- mu[[i]]$summary |> 
    filter(!str_starts(param, "beta_") | str_starts(param, "beta_timepoint"), param != "B_timepoint")
  metamu[[i]]$summary <- metamu[[i]]$summary |> 
    filter(!str_starts(param, "beta_") | str_starts(param, "beta_timepoint"), param != "B_timepoint")
}


check_fits(list(unadjusted$age, mu$age))

devplot(unadjusted$age, mu$age, output="table") |> 
  arrange(-diff)

mu$age$results |> filter(str_starts(param, "mu")) |> slice(72)

mods <- list()

for (i in 1:length(chars)) {
  message(chars[i])
  load_chars(chars[i])
  classes <- distinct(id_lookup, class, cl) |> arrange(cl)
  
  mu[[chars[i]]]$summary <- mu[[chars[i]]]$summary |> filter(!str_starts(param, paste0("beta_mu|beta_", chars[i])), `2.5%` > -8)
  
  mods[[chars[i]]] <- list()
  c <- readRDS(paste0("R/saved-results/metareg/", chars[i], ".rds"))
  
  c$base$summary <- c$base$summary |>
    mutate(
      cl = if_else(
        str_starts(param, "beta_"),
        suppressWarnings(as.numeric(str_sub(param, -2, -2))), NA
      ),
      label = if_else(
        str_starts(param, "beta_"),
        paste0(str_sub(param, 6, -3), classes$class[cl], "]"), label
      )
    ) |> select(-cl) |>
    filter(`97.5%` < 10, !str_starts(param, paste0("beta_", chars[i])))
  c$com$summary <- c$com$summary |>
    mutate(
      cl = if_else(
        str_starts(param, "beta_"),
        suppressWarnings(as.numeric(str_sub(param, -2, -2))), NA
      ),
      label = if_else(
        str_starts(param, "beta_"),
        paste0(str_sub(param, 6, -3), classes$class[cl], "]"), label
      )
    ) |> select(-cl) |>
    filter(`97.5%` < 10, !str_starts(param, paste0("beta_", chars[i])))
  mods[[chars[i]]][[1]] <- c$base
  mods[[chars[i]]][[2]] <- c$com
  mods[[chars[i]]][[3]] <- mu[[chars[i]]]
  names(mods[[chars[i]]]) <- c("unadjusted", "cov", "cov/baseline")
  
  write.csv(compare_jags(mods[[chars[i]]]), paste0("output/metareg/", chars[i], ".csv"))
  forests(mods[[chars[i]]]) + labs(title = chars[i])
  ggsave(paste0("output/metareg/", chars[i], "_mu.png"), height = 10, width = 7)
  
}

lapply(mods, \(x) lapply(x, \(y) y$summary)) |> 
  unlist(recursive = FALSE) |> 
  bind_rows(.id = "model") |> 
  filter(str_starts(param, "B_"), param != "B_mu") |> 
  mutate(model = str_remove(model, ".*?[.]")) |> 
  ggplot(aes(x = mean, y = param)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = model), position = position_dodge(width = 0.5, reverse = TRUE), shape=15) +
  scale_colour_discrete("", labels = c("Covariate adjustment", "Covariate and baseline adjustment")) +
  theme_bw() +
  labs(y = "Parameter", x = "Change in probit response rate") +
  theme(legend.position = "bottom")
  
ggsave("output/metareg/@meta_summary.png", height = 5, width = 5)


