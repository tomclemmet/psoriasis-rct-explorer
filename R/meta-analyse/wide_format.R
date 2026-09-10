library(DBI)
library(dplyr)
library(tidyr)
drug_class_lookup <- read.csv("R/meta-analyse/trt_class.csv")
pasi_drugs <- c(
  "Adalimumab", "Apremilast", "Bimekizumab", "Brodalumab", "Certolizumab", "Cyclosporin",
  "Deucravacitinib", "Etanercept", "Fumaric acid esters", "Guselkumab", "Icotrokinra", "Infliximab",
  "Ixekizumab", "Izokibep", "Methotrexate", "Mirikizumab", "Netakimab", "Placebo",
  "Risankizumab", "Roflumilast", "Secukinumab", "Sonelokimab", "Tildrakizumab", "Tofacitinib",
  "Ustekinumab"
)

gen_pasi_jags <- function(drugs = pasi_drugs, direct = TRUE) {
  
  con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")
  
  pasi_condensed <<- dbReadTable(con, "v_pasi") |> 
    group_by(ref_id) |> 
    ungroup() |> 
    filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |> 
    summarise(.by = c(trial, ref_id, drug, pasi_high_rob, pop_res), n = sum(n), pasi50 = sum(pasi50), 
              pasi75 = sum(pasi75), pasi90 = sum(pasi90), pasi100 = sum(pasi100)) |> 
    filter(pop_res %notin% c("Comorbidity restrictions", "Cormorbidity restrictions", "Systemic-naïve", "Inadequate response to ustekinumab",
                             "Nail psoriasis", "Psoriatic arthritis", "Scalp psoriasis")) |> 
    filter(pasi_high_rob == 0)
  dbDisconnect(con)
  
  if (direct == TRUE) {
    pasi_condensed <<- pasi_condensed |>
      group_by(trial) |> filter(drug %in% drugs) |> filter(n() > 1) |> ungroup()
  } else {
    pasi_condensed <<- pasi_condensed |>
      group_by(trial) |> filter(sum(drug %in% drugs) < length(drugs)) |> ungroup()
  }
  
  if ("Placebo" %in% drugs) {
    pasi_drugs <<- c("Placebo", setdiff(sort(unique(pasi_condensed$drug)), "Placebo"))
  } else {
    pasi_drugs <<- sort(unique(pasi_condensed$drug))
  }
  
  nth_largest <- function(n, ...) {
    vec <- c(...)
    sort(vec, decreasing = TRUE, na.last = TRUE)[n]
  }
  
  nth_non_na <- function(n, ...) {
    vec <- c(...)
    which(!is.na(vec))[n]
  }

  pasi_wide <<- pasi_condensed |> 
    left_join(drug_class_lookup, by = "drug") |> 
    group_by(ref_id) |> 
    mutate(across(pasi50:pasi100, \(x) if (any(is.na(x))) NA else x)) |>
    ungroup() |> 
    mutate(t = as.numeric(factor(drug, levels = pasi_drugs)), 
           cl = as.numeric(factor(class, levels = c("placebo", setdiff(sort(class), "placebo")))), 
           .after = drug) |>
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
    select(-c(drug, class, trial, pasi_high_rob, pop_res, n, pasi50:pasi100)) |> 
    mutate(.by = ref_id, na = n(), arm_no = row_number(t)) |> arrange(ref_id, t) |> 
    pivot_wider(names_from = arm_no, values_from = t:n5, names_glue = "a{arm_no}{.value}") |> 
    relocate(na, .before = nc)
  
  pasi_id_lookup <<- pasi_condensed |> 
    left_join(drug_class_lookup, by = "drug") |> 
    select(trial, ref_id, drug, class, n:pasi100) |>
    mutate(.by = ref_id, arm_no = row_number()) |> 
    filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |> 
    group_by(ref_id) |> 
    mutate(across(pasi50:pasi100, \(x) if (any(is.na(x))) NA else x)) |>
    ungroup() |> 
    mutate(t = as.numeric(factor(drug, levels = pasi_drugs)), 
           cl = as.numeric(factor(class, levels = c("placebo", setdiff(sort(class), "placebo"))))) |> 
    distinct(drug, class, t, cl) |> arrange(t) |> as.data.frame()
  
  row_lookup <<- pasi_wide |> 
    mutate(row = row_number()) |> 
    select(row, ref_id, starts_with("a") & ends_with("t")) |> 
    pivot_longer(starts_with("a") & ends_with("t"), names_to = "arm_no", values_to = "t") |> 
    mutate(arm_no = as.numeric(substr(arm_no, 2, 2))) |> 
    filter(!is.na(t)) |> 
    mutate(drug = pasi_id_lookup$drug[t]) |> 
    left_join(distinct(pasi_condensed, trial, ref_id), by = "ref_id", relationship = "many-to-one")
  
  list(
    ns = nrow(pasi_wide),
    nt = max(pasi_id_lookup$t),
    ncl = max(pasi_id_lookup$cl),
    Cmax = 5,
    mmu = 0.45,
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
      simplify2array()
  )
}

pasi_jags <- gen_pasi_jags()