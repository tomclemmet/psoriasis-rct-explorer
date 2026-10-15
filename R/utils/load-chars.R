# Loads relevant baseline characteristics
library(DBI)
library(dplyr)
library(tidyr)

# Load data
reg_ids <- c()

load_chars <- function(my_chars) {
  
  char_labs <- c()
  for (i in 1:length(my_chars)) {
    char_labs[i] <- switch(
      my_chars[i],
      "age" = "Age",
      "male" = "Sex (n male)",
      "white" = "Ethnicity (n white)",
      "weight" = "Weight",
      "bmi" = "BMI",
      "pasi" = "PASI",
      "duration" = "Duration of psoriasis",
      "psa" = "Psoriatic arthritis",
      "bio" = "Biologic agents"
    )
  }

  con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")
  
  pasi_ids <- dbReadTable(con, "v_pasi") |> 
    filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |>
    # group_by(ref_id) |> filter(any(drug == "Placebo"), any(drug == "Ustekinumab"), drug %in% c("Ustekinumab", "Placebo")) |> ungroup() |>
    pull(ref_id) |> unique()

  chars <- dbReadTable(con, "v_baseline") |> 
    # group_by(ref_id) |> filter(any(drug == "Placebo"), any(drug == "Ustekinumab"), drug %in% c("Ustekinumab", "Placebo")) |> ungroup() |>
    filter(characteristic %in% char_labs, !is.na(k) | !is.na(mean), ref_id %in% pasi_ids) |>
    mutate(arm_no = if_else(!is.na(parent_arm_no), parent_arm_no, arm_no),
           t = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo"))))) |> 
    summarise(
      .by = c(trial, ref_id, drug, t, arm_no, characteristic, data_type),
      sd = sqrt((sum((n - 1) * sd^2) + sum(n * (mean - (sum(mean * n) / sum(n)))^2)) / (sum(n) - 1)),
      mean = sum(mean * n) / sum(n),
      n = sum(n),
      k = sum(k),
    ) |> 
    arrange(ref_id, t, arm_no) |> mutate(.by = ref_id, arm_no = row_number())
    

  dbDisconnect(con)

  reg_ids <<- chars |>
    distinct(ref_id, characteristic) |> 
    group_by(ref_id) |> 
    pull(ref_id) |> unique()

  char_list <- list()

  for (i in 1:length(my_chars)) {
    
    char_list[[my_chars[i]]] <- chars |> 
      filter(ref_id %in% reg_ids, characteristic == char_labs[i]) |> 
      mutate(x = if_else(data_type == "Continuous", mean, k/n)) |> 
      select(ref_id, arm_no, x) |>
      pivot_wider(names_from = arm_no, values_from = x, names_sort = TRUE) |> 
      select(-ref_id) |> 
      as.matrix()
    
  }
  names(char_list) <- my_chars
  source("R/utils/load-pasi-reg.R")
  
  c(jags_inputs, char_list)
}

