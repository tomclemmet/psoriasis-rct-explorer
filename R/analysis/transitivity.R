# Loads relevant baseline characteristics
library(DBI)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(purrr)
library(stringr)
theme_set(theme_bw())

# Load data
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

my_chars <- c(
  "Age", "Sex (n male)", "Weight",  "BMI", "PASI", "Biologic agents", 
    "Ethnicity (n white)", "Duration of psoriasis"
)

pasi_ids <- dbReadTable(con, "v_pasi") |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |>
  pull(ref_id) |> unique()

chars <- dbReadTable(con, "v_baseline") |> 
  arrange(ref_id, arm_no) |> 
  filter(characteristic %in% my_chars, !is.na(k) | !is.na(mean)) |>
  mutate(arm_no = if_else(!is.na(parent_arm_no), parent_arm_no, arm_no)) |> 
  summarise(
    .by = c(trial, ref_id, drug, arm_no, characteristic, data_type),
    sd = sqrt((sum((n - 1) * sd^2) + sum(n * (mean - (sum(mean * n) / sum(n)))^2)) / (sum(n) - 1)),
    mean = sum(mean * n) / sum(n),
    n = sum(n),
    k = sum(k),
  ) |> 
  filter(ref_id %in% pasi_ids) |> 
  mutate(.by = ref_id, arm_no = dense_rank(arm_no))

dbDisconnect(con)

chars |> 
  filter(characteristic == "Age") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Mean age (yrs)")
# Narrow range and mostly overlapping, no concerns

chars |> 
  filter(characteristic == "Sex (n male)") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Sex (% male)")
# Consistent rate of around 70%, no concerns

chars |> 
  filter(characteristic == "Weight") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot(aes(alpha = n)) + 
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16)
  labs(x = "Mean weight (kg)")
# Signs of heterogeneity but no major transitivity concerns

chars |> 
  filter(characteristic == "BMI") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Mean BMI (kg/m^2)")
# Similar to weight

chars |> 
  filter(characteristic == "PASI") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Mean PASI")
# Lots of outliers but main comparisons are consistent

chars |> 
  filter(characteristic == "Biologic agents") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair)), 
         yr = str_sub(trial, -4, -1), yr = as.numeric(if_else(yr == "004a", "2004", if_else(yr=="4309", "2022", yr)))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id, colour = yr), alpha = 0.7, shape = 16) +
  scale_color_viridis_c() +
  labs(x = "Biologic experienced (%)")
# Biggest transitivity concern, wide range and differences by comparison

chars |> 
  filter(characteristic == "Ethnicity (n white)") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Ethnicity (% white)")
# Narrow range, no major concerns

chars |> 
  filter(characteristic == "Duration of psoriasis") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot(aes(alpha = n)) +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Duration of psoriais (yrs)")
# Narrow range, no major concern
  
# PLAN: Adjust for weight and biologics