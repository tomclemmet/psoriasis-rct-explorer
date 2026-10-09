library(DBI)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

# Extract PASI response data in long format
pasi_long <- dbReadTable(con, "v_pasi") |> 
  select(trial, ref_id, drug, arm_no, n:pasi100) |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |>
  pivot_longer(pasi50:pasi100, names_to = "code", values_to = "k") |> 
  filter(!is.na(k)) |> 
  mutate(subgroup = "0verall")

subgrps <- dbReadTable(con, "v_subgroups") |> 
  select(trial, ref_id, arm_no, drug, subgroup, code, k, n) |> 
  filter(code %in% c("pasi50", "pasi75", "pasi90", "pasi100"))

comb <- pasi_long |> 
  filter(ref_id %in% unique(subgrps$ref_id)) |> 
  bind_rows(subgrps) |> 
  mutate(p = k/n, code = factor(code, levels = c("pasi50", "pasi75", "pasi90", "pasi100")),
         subgroup = if_else(subgroup %in% c("Ever biologic"), "Biologic-experienced", subgroup),
         subgroup = if_else(subgroup %in% c("Never biologic"), "Biologic-naive", subgroup)) |>
  mutate(lab = paste0("Arm ", dense_rank(arm_no), ": ", drug), 
         .by = c(ref_id, subgroup)) |> 
  group_by(ref_id, code) |> 
  filter(any(subgroup %in% c("Biologic-naive", "Biologic-experienced"))) |> 
  ungroup() |> 
  filter(subgroup %in% c("0verall", "Biologic-naive", "Biologic-experienced"))

comb |> 
  filter(ref_id %in% unique(comb$ref_id)[5]) |> 
  ggplot() +
  geom_col(aes(x = p, y = subgroup, fill = code), position = "identity") +
  facet_wrap(~ lab, scale = "free_x")

comb |> 
  filter(subgroup != "0verall", code == "pasi75") |> 
  select(-k, -n) |> 
  pivot_wider(names_from = subgroup, values_from = p) |> 
  ggplot(aes(x = `Biologic-experienced`, y = `Biologic-naive`)) +
  geom_point(aes(colour = drug)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_viridis_d(option = "turbo")
ggsave("output/subgroup_bio.png", height = 6, width = 8)

