library(DBI)
library(dplyr)
library(tidyr)
class_lookup <- read.csv("R/meta-analyse/trt_class.csv")

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

pasi <- dbReadTable(con, "v_pasi")
dbDisconnect(con)

drug_order <- c(
  "Placebo", "Acitretin", "Adalimumab", "Apremilast", "Bimekizumab",
  "Brodalumab", "Certolizumab", "Cyclosporin", "Deucravacitinib", "Etanercept",
  "Fumaric acid esters", "Guselkumab", "Icotrokinra", "Infliximab", "Ixekizumab",
  "Izokibep", "Methotrexate", "Mirikizumab", "Netakimab", "Orismilast",
  "Phototherapy", 
  "Risankizumab", "Roflumilast", "Secukinumab", "Sonelokimab",
  "Tildrakizumab", "Tofacitinib", "Ustekinumab", "Xeligekimab", "Zasocitinib"
)

nth_largest <- function(n, ...) {
  vec <- c(...)
  sort(vec, decreasing = TRUE, na.last = TRUE)[n]
}

nth_non_na <- function(n, ...) {
  vec <- c(...)
  which(!is.na(vec))[n]
}

pasi_wide <- pasi |> 
  select(trial, ref_id, arm_no, drug, n:pasi100) |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |> 
  # group_by(trial, ref_id, drug) |> 
  # summarise(n = sum(n), pasi50 = sum(pasi50), pasi75 = sum(pasi75), 
  #           pasi90 = sum(pasi90), pasi100 = sum(pasi100), .groups = "drop") |> 
  group_by(ref_id) |> 
  mutate(across(pasi50:pasi100, \(x) if (any(is.na(x))) NA else x)) |>
  ungroup() |> 
  mutate(t = as.numeric(factor(drug, levels = drug_order)), .after = drug) |>
  rowwise() |> mutate(
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
    C1 = 1,
    C2 = nth_non_na(1, pasi50, pasi75, pasi90, pasi100) + 1,
    C3 = nth_non_na(2, pasi50, pasi75, pasi90, pasi100) + 1,
    C4 = nth_non_na(3, pasi50, pasi75, pasi90, pasi100) + 1,
    C5 = nth_non_na(4, pasi50, pasi75, pasi90, pasi100) + 1,
    nc = sum(!is.na(c(C1, C2, C3, C4, C5)))
  ) |> ungroup() |> relocate(C1:nc, .before = t) |> 
  select(-c(drug, trial, n, pasi50:pasi100)) |> 
  mutate(.by = ref_id, na = n(), arm_no = row_number(t)) |> arrange(ref_id, t) |> 
  pivot_wider(names_from = arm_no, values_from = t:n5, names_glue = "a{arm_no}{.value}") |> 
  relocate(na, .before = nc)

pasi_jags <- list(
  ns = nrow(pasi_wide),
  nt = n_distinct(filter(pasi, !if_all(pasi50:pasi100, \(x) is.na(x)))$drug),
  Cmax = max(pasi_wide$nc),
  mmu = 0.6,
  na = pasi_wide$na,
  nc = pasi_wide$nc,
  t = select(pasi_wide, a1t:a5t) |> as.matrix(),
  C = select(pasi_wide, C1:C5) |> as.matrix(),
  r = list(as.matrix(select(pasi_wide, a1r1:a5r1)),
           as.matrix(select(pasi_wide, a1r2:a5r2)),
           as.matrix(select(pasi_wide, a1r3:a5r3)),
           as.matrix(select(pasi_wide, a1r4:a5r4)),
           as.matrix(select(pasi_wide, a1r5:a5r5))) |> 
    simplify2array(),
  n = list(as.matrix(select(pasi_wide, a1n1:a5n1)),
           as.matrix(select(pasi_wide, a1n2:a5n2)),
           as.matrix(select(pasi_wide, a1n3:a5n3)),
           as.matrix(select(pasi_wide, a1n4:a5n4)),
           as.matrix(select(pasi_wide, a1n5:a5n5))) |> 
    simplify2array()
)


# source("R/meta-analyse/pasi-jags-nma.R")
# source("R/meta-analyse/ma-utils.R")
# m <- pso_jags(pasi_jags, filename = "JAGS/re_fez_u.jags", effects = "random", cutpoints = "fixed", baseline = "unadjusted") 
# process_jags(m)$summary
