source("R/utils/plots.R")
source("R/utils/load-chars.R")
classes <- distinct(id_lookup, class, cl) |> arrange(cl)
char <- "bio"

m <- readRDS(paste0("R/saved-results/metareg/", char, ".rds"))
load_chars(char)

check_fits(m)

# for (i in 1:length(m)) {
#   message(names(m)[i])
#   f <- m[[i]]$summary |> 
#     filter(Rhat > 1.1)
#   if(nrow(f)>0) print(f)
# }

for (i in 1:length(m)) {
  m[[i]]$summary <- m[[i]]$summary |> 
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
    filter(`97.5%` < 10)
}

write.csv(compare_jags(m), paste0("output/metareg/", char, ".csv"))
forests(m)
ggsave(paste0("output/metareg/", char, ".png"), height = 10, width = 7)
