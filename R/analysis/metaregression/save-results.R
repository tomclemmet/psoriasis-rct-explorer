source("R/utils/plots.R")
source("R/utils/load-chars.R")
char <- "bio_mu"

m <- readRDS(paste0("R/saved-results/metareg/", char, ".rds"))
load_chars("bio")
# source("R/utils/load-pasi.R")
classes <- distinct(id_lookup, class, cl) |> arrange(cl)

check_fits(m)

for (i in 1:length(m)) {
  message(names(m)[i])
  f <- m[[i]]$summary |>
    filter(Rhat > 1.1)
  if(nrow(f)>0) print(f)
}

# for (i in 1:length(m)) {
#   m[[i]]$summary <- m[[i]]$summary |> 
#     mutate(
#       cl = if_else(
#         str_starts(param, "beta_"),
#         suppressWarnings(as.numeric(str_sub(param, -2, -2))), NA
#       ),
#       label = if_else(
#         str_starts(param, "beta_"), 
#         paste0(str_sub(param, 6, -3), classes$class[cl], "]"), label
#       )
#     ) |> select(-cl) |> 
#     filter(`97.5%` < 10)
# }
for (i in 1:length(m)) {
  m[[i]]$summary <- m[[i]]$summary |> 
    filter(!str_starts(param, "beta_mu"))
}
m$com$summary <- m$com$summary |> 
  mutate(mean = if_else(label != "timepoint[tsm]" | is.na(label), mean, 0))


write.csv(compare_jags(m), paste0("output/metareg/", char, ".csv"))
forests(m)
ggsave(paste0("output/metareg/", char, ".png"), height = 10, width = 7)
