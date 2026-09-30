library(DBI)
library(dplyr)
library(tidyr)
drug_class_lookup <- read.csv("R/trt_class.csv")

# Load data
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

# List drugs for analysis (allows filtering)
pasi_drugs <- dbReadTable(con, "v_pasi") |> 
  distinct(drug) |> 
  # filter(drug != "Izokibep") |> # CHANGE WHEN PhIII trial included
  mutate(drug = c("Placebo", sort(setdiff(drug, "Placebo")))) |> 
  pull(drug)

# Extract PASI response data in long format
pasi_long <- dbReadTable(con, "v_pasi") |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |>
  mutate(arm_no = if_else(!is.na(parent_arm_no), parent_arm_no, arm_no)) |> 
  summarise(
    .by = c(trial, ref_id, drug, pop_res, pasi_high_rob, timepoint, timepoint_unit),
    n = sum(n), pasi50 = sum(pasi50), pasi75 = sum(pasi75), 
    pasi90 = sum(pasi90), pasi100 = sum(pasi100)
  ) |> 
  mutate(timepoint = if_else(timepoint == 4 & timepoint_unit == "mo", 16, timepoint),
         pop_res = na_if(pop_res, "null")) |>
  select(-timepoint_unit)
  # summarise(.by = c(trial, ref_id, drug, pasi_high_rob, pop_res, timepoint), n = sum(n), pasi50 = sum(pasi50), 
            # pasi75 = sum(pasi75), pasi90 = sum(pasi90), pasi100 = sum(pasi100))
  # filter(pop_res %notin% c("Inadequate response to ustekinumab")) |> 
  # filter(pasi_high_rob == 0)

dbDisconnect(con)

# List timepoints
pasi_long |> 
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
pasi_wide <- pasi_long |> 
  left_join(drug_class_lookup, by = "drug") |> # Add drug classes
  group_by(ref_id) |> # Align trials where the arms report results for different cutpoints
  mutate(across(pasi50:pasi100, \(x) if (any(is.na(x))) NA else x)) |>
  ungroup() |> 
  mutate(t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))), 
         cl = as.numeric(factor(class, levels = c("placebo", setdiff(sort(class), "placebo")))), 
         .after = drug) |> # Assign treatment and class ids
  rowwise() |> mutate( # Calculate conditional counts
    r1 = n - nth_largest(1, pasi50, pasi75, pasi90, pasi100),
    r2 = nth_largest(1, pasi50, pasi75, pasi90, pasi100) - 
      coalesce(nth_largest(2, pasi50, pasi75, pasi90, pasi100), 0),
    r3 = nth_largest(2, pasi50, pasi75, pasi90, pasi100) - 
      coalesce(nth_largest(3, pasi50, pasi75, pasi90, pasi100), 0),
    r4 = nth_largest(3, pasi50, pasi75, pasi90, pasi100) - 
      coalesce(nth_largest(4, pasi50, pasi75, pasi90, pasi100), 0),
    r5 = nth_largest(4, pasi50, pasi75, pasi90, pasi100), 
    n1 = n,
    n2 = if_else(is.na(r2), NA, n1 - coalesce(r1, 0)),
    n3 = if_else(is.na(r3), NA, n2 - coalesce(r2, 0)),
    n4 = if_else(is.na(r4), NA, n3 - coalesce(r3, 0)),
    n5 = if_else(is.na(r5), NA, n4 - coalesce(r4, 0)),
    C1 = 1, # Record categories
    C2 = nth_non_na(1, pasi50, pasi75, pasi90, pasi100) + 1,
    C3 = nth_non_na(2, pasi50, pasi75, pasi90, pasi100) + 1,
    C4 = nth_non_na(3, pasi50, pasi75, pasi90, pasi100) + 1,
    C5 = nth_non_na(4, pasi50, pasi75, pasi90, pasi100) + 1,
    nc = sum(!is.na(c(C1, C2, C3, C4, C5)))
  ) |> ungroup() |> relocate(C1:nc, .before = t) |> # Move/drop non-pivot columns
  select(-c(drug, class, trial, pasi_high_rob, pop_res, timepoint, n, pasi50:pasi100)) |> 
  mutate(.by = ref_id, na = n(), arm_no = dense_rank(t)) |> arrange(ref_id, arm_no) |> # Add arm info
  pivot_wider(names_from = arm_no, values_from = t:n5, names_glue = "a{arm_no}{.value}") |> # Pivot to wide format
  relocate(na, .before = nc)

# Treatment and class IDs
pasi_id_lookup <- pasi_long |> 
  left_join(drug_class_lookup, by = "drug") |> 
  select(trial, ref_id, drug, class, n:pasi100) |>
  mutate(.by = ref_id, arm_no = row_number()) |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |> 
  group_by(ref_id) |> 
  mutate(across(pasi50:pasi100, \(x) if (any(is.na(x))) NA else x)) |>
  ungroup() |> 
  mutate(t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))), 
         cl = as.numeric(factor(class, levels = c("placebo", setdiff(sort(class), "placebo"))))) |> 
  distinct(drug, class, t, cl) |> arrange(t) |> as.data.frame()

# IDs by pasi_wide row number
row_lookup <- pasi_wide |> 
  mutate(row = row_number()) |> 
  select(row, ref_id, starts_with("a") & ends_with("t")) |> 
  pivot_longer(starts_with("a") & ends_with("t"), names_to = "arm_no", values_to = "t") |> 
  mutate(arm_no = as.numeric(substr(arm_no, 2, 2))) |> 
  filter(!is.na(t)) |> 
  mutate(drug = pasi_id_lookup$drug[t]) |> 
  left_join(distinct(pasi_long, trial, ref_id), by = "ref_id", relationship = "many-to-one")

# Collect data in JAGS format
pasi_jags <- list(
  ns = nrow(pasi_wide),
  nt = max(pasi_id_lookup$t),
  ncl = max(pasi_id_lookup$cl),
  Cmax = 5,
  na = pasi_wide$na,
  nc = pasi_wide$nc,
  t = select(pasi_wide, starts_with("a") & ends_with("t")) |> as.matrix(),
  cl = pasi_id_lookup$cl,
  C = select(pasi_wide, C1:C5) |> as.matrix(),
  r = list(as.matrix(select(pasi_wide, ends_with("r1"))),
           as.matrix(select(pasi_wide, ends_with("r2"))),
           as.matrix(select(pasi_wide, ends_with("r3"))),
           as.matrix(select(pasi_wide, ends_with("r4"))),
           as.matrix(select(pasi_wide, ends_with("r5")))) |> 
    simplify2array(),
  n = list(as.matrix(select(pasi_wide, ends_with("n1"))),
           as.matrix(select(pasi_wide, ends_with("n2"))),
           as.matrix(select(pasi_wide, ends_with("n3"))),
           as.matrix(select(pasi_wide, ends_with("n4"))),
           as.matrix(select(pasi_wide, ends_with("n5")))) |> 
    simplify2array(),
  timepoint = distinct(pasi_long, ref_id, timepoint)$timepoint
)

