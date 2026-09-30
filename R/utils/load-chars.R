# Loads relevant baseline characteristics
library(DBI)
library(dplyr)
library(tidyr)

# Load data
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

my_chars <- c(
  "Weight", "Biologic agents"
)

pasi_ids <- dbReadTable(con, "v_pasi") |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |>
  pull(ref_id) |> unique()

chars <- dbReadTable(con, "v_baseline") |> 
  mutate(arm_no = if_else(!is.na(parent_arm_no), parent_arm_no, arm_no),
         t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo"))))) |> 
  mutate(.by = ref_id, arm_no = dense_rank(t)) |> 
  arrange(ref_id, characteristic, arm_no) |> 
  filter(characteristic %in% my_chars, !is.na(k) | !is.na(mean)) |>
  summarise(
    .by = c(trial, ref_id, drug, arm_no, characteristic, data_type),
    sd = sqrt((sum((n - 1) * sd^2) + sum(n * (mean - (sum(mean * n) / sum(n)))^2)) / (sum(n) - 1)),
    mean = sum(mean * n) / sum(n),
    n = sum(n),
    k = sum(k),
  ) |> 
  filter(ref_id %in% pasi_ids)

dbDisconnect(con)

reg_ids <- chars |>
  distinct(ref_id, characteristic) |> 
  group_by(ref_id) |> 
  filter(
    # any(characteristic == "Sex (n male)"),
    any(characteristic == "Biologic agents"),
    any(characteristic == "Weight"),
    # any(characteristic == "Age"),
    # any(characteristic == "PASI"),
    # any(characteristic == "Ethnicity (n white)"),
    # any(characteristic == "Duration of psoriasis")
  ) |> 
  pull(ref_id) |> unique()

char_list <- list()

for (i in 1:length(my_chars)) {
  
  char_list[[my_chars[i]]] <- chars |> 
    filter(ref_id %in% reg_ids, characteristic == my_chars[i]) |> 
    mutate(x = if_else(data_type == "Continuous", mean, k/n)) |> 
    select(ref_id, arm_no, x) |>
    pivot_wider(names_from = arm_no, values_from = x, names_sort = TRUE) |> 
    select(-ref_id) |> 
    as.matrix()
  
}
names(char_list) <- c("weight", "bio")

char_list$weight <- char_list$weight / 10

