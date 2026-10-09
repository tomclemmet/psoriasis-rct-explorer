# Loads relevant baseline characteristics
library(DBI)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(purrr)
library(stringr)
theme_set(theme_bw())
source("R/utils/load-chars.R")

# Load data
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

my_chars <- c(
  "Age", "Sex (n male)", "Weight",  "BMI", "PASI", "Biologic agents", 
    "Ethnicity (n white)", "Duration of psoriasis", "Psoriatic arthritis"
)

pasi_ids <- dbReadTable(con, "v_pasi") |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |>
  pull(ref_id) |> unique()

chars <- dbReadTable(con, "v_baseline") |> 
  filter(characteristic %in% my_chars, !is.na(k) | !is.na(mean), ref_id %in% pasi_ids) |>
  mutate(arm_no = if_else(!is.na(parent_arm_no), parent_arm_no, arm_no),
         t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo"))))) |> 
  summarise(
    .by = c(trial, ref_id, drug, t, arm_no, characteristic, data_type),
    sd = sqrt((sum((n - 1) * sd^2) + sum(n * (mean - (sum(mean * n) / sum(n)))^2)) / (sum(n) - 1)),
    mean = sum(mean * n) / sum(n),
    n = sum(n),
    k = sum(k),
  ) |> 
  arrange(ref_id, t, arm_no) |> mutate(.by = c(ref_id, characteristic), arm_no = row_number())

dbDisconnect(con)

chars |> 
  filter(characteristic == "Age") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Mean age (yrs)")
# Narrow range and mostly overlapping, no concerns
ggsave("output/transitivity/age.png", height = 10, width = 8)

chars |> 
  filter(characteristic == "Sex (n male)") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Sex (% male)")
# Consistent rate of around 70%, no concerns
ggsave("output/transitivity/sex.png", height = 10, width = 8)


chars |> 
  filter(characteristic == "Weight") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot() + 
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16)
  labs(x = "Mean weight (kg)")
# Signs of heterogeneity but no major transitivity concerns
ggsave("output/transitivity/weight.png", height = 10, width = 8)

chars |> 
  filter(characteristic == "BMI") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Mean BMI (kg/m^2)")
# Similar to weight
ggsave("output/transitivity/bmi.png", height = 10, width = 8)


chars |> 
  filter(characteristic == "PASI") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Mean PASI")
# Lots of outliers but main comparisons are consistent
ggsave("output/transitivity/pasi.png", height = 10, width = 8)

chars |> 
  filter(characteristic == "Psoriatic arthritis") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "With psoriatic arthritis (%)")
# Some concerns but mostly in single-trial comparisons
ggsave("output/transitivity/psa.png", height = 10, width = 8)

chars |> 
  filter(characteristic == "Ethnicity (n white)") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Ethnicity (% white)")
# Narrow range, no major concerns
ggsave("output/transitivity/ethnicity.png", height = 10, width = 8)


chars |> 
  filter(characteristic == "Duration of psoriasis") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair))) |> 
  ggplot(aes(x = mean, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id), alpha = 0.5, shape = 16) +
  labs(x = "Duration of psoriais (yrs)")
# Narrow range, no major concern
ggsave("output/transitivity/duration.png", height = 10, width = 8)

chars |> 
  filter(characteristic == "Biologic agents") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair)), 
         yr = str_sub(trial, -4, -1), yr = as.numeric(if_else(yr == "004a", "2004", if_else(yr=="4309", "2022", yr)))) |> 
  ggplot(aes(x = k/n, y = pair)) +
  geom_boxplot() +
  # geom_point(aes(size = n, text = ref_id, colour = yr), alpha = 0.7, shape = 16) +
  # scale_color_viridis_c() +
  labs(x = "Biologic experienced (%)")
# Biggest transitivity concern, wide range and differences by comparison
ggsave("output/transitivity/bio.png", height = 10, width = 8)


# Add PASI to graph
load_chars("bio")
ggplotly(chars |> 
  filter(characteristic == "Biologic agents") |> 
  mutate(comps = list(unique(drug)), .by = ref_id) |> 
  mutate(pair = map(comps, \(x) combn(x, 2, simplify = FALSE))) |> 
  unnest(pair) |> 
  mutate(na = n_distinct(ref_id),
         pair = paste0(map_chr(pair, paste, collapse = " vs "), " (", na, " trials)"),
         .by = pair) |> 
  filter(str_detect(pair, drug)) |> 
  arrange(na, pair) |> 
  mutate(pair = factor(pair, levels = unique(pair)), 
         yr = str_sub(trial, -4, -1), yr = as.numeric(if_else(yr == "004a", "2004", if_else(yr=="4309", "2022", yr)))) |> 
  left_join(select(pasi_long, ref_id, arm_no, pasi75), by = c("ref_id", "arm_no")) |> 
  summarise(.by = c(pair, drug, ref_id), pasi75 = sum(pasi75, na.rm = TRUE), n = sum(n), k = sum(k)) |> 
  mutate(.by = c(pair, ref_id), t = paste0("t", row_number())) |> select(-drug) |> 
  pivot_wider(names_from = t, values_from = k:pasi75) |> 
  mutate(rd = pasi75_t2 / n_t2 - pasi75_t1 / n_t1,
         p = (k_t1 + k_t2) / (n_t1 + n_t2), .keep = "unused") |>
  filter(rd != 0) |> 
  ggplot(aes(x = rd, y = pair)) +
  geom_point(aes(colour = p, text = ref_id)) +
  scale_color_viridis_c() +
  labs(x = "Biologic experienced (%)"), tooltip = "text")

chars |> 
  filter(characteristic == "Biologic agents") |> 
  group_by(ref_id) |> 
  filter(any(drug == "Ustekinumab"), any(drug == "Placebo"), drug %in% c("Ustekinumab", "Placebo")) |> 
  ungroup() |> 
  left_join(select(pasi_long, ref_id, arm_no, pasi75), by = c("ref_id", "arm_no")) |> 
  summarise(.by = c(drug, ref_id), pasi75 = sum(pasi75, na.rm = TRUE), n = sum(n), k = sum(k)) |> 
  arrange(ref_id, desc(drug)) |> 
  mutate(.by = c(ref_id), t = paste0("t", row_number())) |> select(-drug) |> 
  pivot_wider(names_from = t, values_from = pasi75:k) |> 
  mutate(rd = pasi75_t1 / n_t1 - pasi75_t2 / n_t2,
         p = (k_t1 + k_t2) / (n_t1 + n_t2), .keep = "unused", 
         n = n_t1 + n_t2) |>
  ggplot(aes(x = p, y = rd)) +
  geom_point(aes(size = n), show.legend = FALSE) +
  geom_smooth(aes(weight = n), method = "lm", se = FALSE) +
  labs(title = "Ustekinumab vs Placebo", x = "Biologic experienced (%)", y = "PASI 75 risk difference")
ggsave("output/bioplot/ust.png", height = 4, width = 4)

chars |> 
  filter(characteristic == "Biologic agents") |> 
  group_by(ref_id) |> 
  filter(any(drug == "Etanercept"), any(drug == "Placebo"), drug %in% c("Etanercept", "Placebo")) |> 
  ungroup() |> 
  left_join(select(pasi_long, ref_id, arm_no, pasi75), by = c("ref_id", "arm_no")) |> 
  summarise(.by = c(drug, ref_id), pasi75 = sum(pasi75, na.rm = TRUE), n = sum(n), k = sum(k)) |> 
  arrange(ref_id, desc(drug)) |> 
  mutate(.by = c(ref_id), t = paste0("t", row_number())) |> select(-drug) |> 
  pivot_wider(names_from = t, values_from = pasi75:k) |> 
  mutate(rd = pasi75_t1 / n_t1 - pasi75_t2 / n_t2,
         p = (k_t1 + k_t2) / (n_t1 + n_t2), .keep = "unused", 
         n = n_t1 + n_t2) |>
  ggplot(aes(x = p, y = rd)) +
  geom_point(aes(size = n), show.legend = FALSE) +
  geom_smooth(aes(weight = n), method = "lm", se = FALSE) +
  labs(title = "Etanercept vs Placebo", x = "Biologic experienced (%)", y = "PASI 75 risk difference")
ggsave("output/bioplot/etn.png", height = 4, width = 4)

chars |> 
  filter(characteristic == "Biologic agents") |> 
  group_by(ref_id) |> 
  filter(any(drug == "Secukinumab"), any(drug == "Placebo"), drug %in% c("Secukinumab", "Placebo")) |> 
  ungroup() |> 
  left_join(select(pasi_long, ref_id, arm_no, pasi75), by = c("ref_id", "arm_no")) |> 
  summarise(.by = c(drug, ref_id), pasi75 = sum(pasi75, na.rm = TRUE), n = sum(n), k = sum(k)) |> 
  arrange(ref_id, desc(drug)) |> 
  mutate(.by = c(ref_id), t = paste0("t", row_number())) |> select(-drug) |> 
  pivot_wider(names_from = t, values_from = pasi75:k) |> 
  mutate(rd = pasi75_t1 / n_t1 - pasi75_t2 / n_t2,
         p = (k_t1 + k_t2) / (n_t1 + n_t2), .keep = "unused", 
         n = n_t1 + n_t2) |>
  filter(rd != 0) |> 
  ggplot(aes(x = p, y = rd)) +
  geom_point(aes(size = n), show.legend = FALSE) +
  geom_smooth(aes(weight = n), method = "lm", se = FALSE) +
  labs(title = "Secukinumab vs Placebo", x = "Biologic experienced (%)", y = "PASI 75 risk difference")  
ggsave("output/bioplot/sec.png", height = 4, width = 4)

chars |> 
  filter(characteristic == "Biologic agents") |> 
  group_by(ref_id) |> 
  filter(any(drug == "Risankizumab"), any(drug == "Placebo"), drug %in% c("Risankizumab", "Placebo")) |> 
  ungroup() |> 
  left_join(select(pasi_long, ref_id, arm_no, pasi75), by = c("ref_id", "arm_no")) |> 
  summarise(.by = c(drug, ref_id), pasi75 = sum(pasi75, na.rm = TRUE), n = sum(n), k = sum(k)) |> 
  arrange(ref_id, desc(drug)) |> 
  mutate(.by = c(ref_id), t = paste0("t", row_number())) |> select(-drug) |> 
  pivot_wider(names_from = t, values_from = pasi75:k) |> 
  mutate(rd = pasi75_t1 / n_t1 - pasi75_t2 / n_t2,
         p = (k_t1 + k_t2) / (n_t1 + n_t2), .keep = "unused", 
         n = n_t1 + n_t2) |>
  filter(rd != 0) |> 
  ggplot(aes(x = p, y = rd)) +
  geom_point(aes(size = n), show.legend = FALSE) +
  geom_smooth(aes(weight = n), method = "lm", se = FALSE) +
  labs(title = "Risankizumab vs Placebo", x = "Biologic experienced (%)", y = "PASI 75 risk difference")  
ggsave("output/bioplot/ris.png", height = 4, width = 4)
