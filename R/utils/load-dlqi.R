library(DBI)
library(dplyr)
library(tidyr)
drug_class_lookup <- read.csv("R/trt_class.csv")

# Load data
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

# List drugs for analysis (allows filtering)
dlqi_drugs <- dbReadTable(con, "v_dlqi") |> 
  distinct(drug) |> 
  # filter(drug != "Izokibep") |> # CHANGE WHEN PhIII trial included
  mutate(drug = c("Placebo", sort(setdiff(drug, "Placebo")))) |> 
  pull(drug)

# Extract DLQI response data in long format
dlqi_long <- dbReadTable(con, "v_dlqi") |> as_tibble() |> 
  filter(!if_all(dlqi_0_1:dlqi_0, \(x) is.na(x))) |>
  mutate(arm_no = if_else(!is.na(parent_arm_no), parent_arm_no, arm_no),
         t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo"))))) |> 
  summarise(
    .by = c(trial, ref_id, drug, t, arm_no, timepoint, timepoint_unit),
    n = sum(n), dlqi_0_1 = sum(dlqi_0_1), dlqi_0 = sum(dlqi_0)
  ) |> 
  arrange(ref_id, t, arm_no) |> mutate(.by = ref_id, arm_no = row_number()) |> 
  mutate(timepoint = if_else(timepoint == 4 & timepoint_unit == "mo", 16, timepoint)) |> 
  select(-timepoint_unit)

dbDisconnect(con)

# List timepoints
dlqi_long |> 
  summarise(
    .by = timepoint, 
    n_pat = sum(n),
    n_arm = n(),
    n_trial = n_distinct(ref_id),
    n_drug = n_distinct(drug)
  ) |> 
  arrange(timepoint)

# Helper function returning the nth largest item from a vector
nth_largest <- function(n, ...) {
  vec <- c(...)
  sort(vec, decreasing = TRUE, na.last = TRUE)[n]
}

# Helper function returning the nth non-na item in a vector
nth_non_na <- function(n, ...) {
  vec <- c(...)
  which(!is.na(vec))[n]
}

# Produce wide-format data frame suitable for JAGS models
dlqi_wide <- dlqi_long |> 
  mutate(.by = ref_id, across(dlqi_0_1:dlqi_0, \(x) if (any(is.na(x))) NA else x)) |> # Align trials where the arms report results for different cutpoints
  rowwise() |> mutate( # Calculate conditional counts
    r1 = n - nth_largest(1, dlqi_0_1, dlqi_0),
    r2 = nth_largest(1, dlqi_0_1, dlqi_0) - 
      coalesce(nth_largest(2, dlqi_0_1, dlqi_0), 0),
    r3 = nth_largest(2, dlqi_0_1, dlqi_0),
    n1 = n,
    n2 = if_else(is.na(r2), NA, n1 - coalesce(r1, 0)),
    n3 = if_else(is.na(r3), NA, n2 - coalesce(r2, 0)),
    C1 = 1, # Record categories
    C2 = nth_non_na(1, dlqi_0_1, dlqi_0) + 1,
    C3 = nth_non_na(2, dlqi_0_1, dlqi_0) + 1,
    nc = sum(!is.na(c(C1, C2, C3)))
  ) |> ungroup() |> relocate(C1:nc, arm_no, .before = t) |> # Move non-pivot columns to before t
  select(-c(drug, trial, timepoint, n, dlqi_0_1:dlqi_0)) |> # Drop additional columns
  mutate(.by = ref_id, na = n()) |> # Add number of arms
  arrange(ref_id, arm_no) |> # Arrange by trial and arm
  pivot_wider(names_from = arm_no, values_from = t:n3, names_glue = "a{arm_no}{.value}") |> # Pivot to wide format
  relocate(na, .before = nc)

# Treatment and class IDs
id_lookup <- dlqi_long |> 
  left_join(drug_class_lookup, by = "drug") |> 
  select(trial, ref_id, drug, class, n:dlqi_0) |>
  mutate(.by = ref_id, arm_no = row_number()) |> 
  filter(!if_all(dlqi_0_1:dlqi_0, \(x) is.na(x))) |> 
  group_by(ref_id) |> 
  mutate(across(dlqi_0_1:dlqi_0, \(x) if (any(is.na(x))) NA else x)) |>
  ungroup() |> 
  mutate(t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))), 
         cl = as.numeric(factor(class, levels = c("placebo", setdiff(sort(class), "placebo"))))) |> 
  distinct(drug, class, t, cl) |> arrange(t) |> as.data.frame()

# IDs by dlqi_wide row number
row_lookup <- dlqi_wide |> 
  mutate(row = row_number()) |> 
  select(row, ref_id, starts_with("a") & ends_with("t")) |> 
  pivot_longer(starts_with("a") & ends_with("t"), names_to = "arm_no", values_to = "t") |> 
  mutate(arm_no = as.numeric(substr(arm_no, 2, 2))) |> 
  filter(!is.na(t)) |> 
  mutate(drug = id_lookup$drug[t]) |> 
  left_join(distinct(dlqi_long, trial, ref_id), by = "ref_id", relationship = "many-to-one")

# Collect data in JAGS format
jags_inputs <- list(
  ns = nrow(dlqi_wide),
  nt = max(id_lookup$t),
  ncl = max(id_lookup$cl),
  Cmax = 3,
  na = dlqi_wide$na,
  nc = dlqi_wide$nc,
  t = select(dlqi_wide, starts_with("a") & ends_with("t")) |> as.matrix(),
  cl = id_lookup$cl,
  C = select(dlqi_wide, C1:C3) |> as.matrix(),
  r = list(as.matrix(select(dlqi_wide, ends_with("r1"))),
           as.matrix(select(dlqi_wide, ends_with("r2"))),
           as.matrix(select(dlqi_wide, ends_with("r3")))) |> 
    simplify2array(),
  n = list(as.matrix(select(dlqi_wide, ends_with("n1"))),
           as.matrix(select(dlqi_wide, ends_with("n2"))),
           as.matrix(select(dlqi_wide, ends_with("n3")))) |> 
    simplify2array(),
  timepoint = distinct(dlqi_long, ref_id, timepoint)$timepoint
)

