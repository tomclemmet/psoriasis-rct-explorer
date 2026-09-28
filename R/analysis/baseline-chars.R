rm(list = ls())
library(dplyr)
library(DBI)
library(ggplot2)
theme_set(theme_bw())

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")
chars <- dbReadTable(con, "v_baseline")

outcomes <- c(
  "Age", "Sex (n male)", "Ethnicity (n white)", "Weight", "BMI", 
  "Duration of psoriasis", "PASI", "Psoriatic arthritis", 
  "Systemic agent", "Biologic agents"
)

chars |> 
  select(ref_id, trial, drug, dose, characteristic, data_type, k, n, mean, sd) |> 
  mutate(mean = if_else(is.na(mean) & !is.na(k), k / n, mean),
         data_type = if_else(characteristic %in% c("PASI", "Duration of psoriasis"), "Continuous-2", data_type)) |> 
  filter(characteristic %in% c(
    "Age", "Sex (n male)", "Ethnicity (n white)", "Weight", "BMI", 
    "Duration of psoriasis", "PASI", "Psoriatic arthritis", 
    "Systemic agent", "Biologic agents"
  )) |> 
  ggplot() +
  geom_boxplot(aes(x = mean, y = characteristic, fill = characteristic), show.legend = FALSE) +
  facet_wrap(~ data_type, scales = "free", space = "free_y", nrow = 3)
ggsave("output/chars.png", height = 7, width = 7)


