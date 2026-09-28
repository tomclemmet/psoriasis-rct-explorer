# Loads relevant baselin characteristics



con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")
chars <- dbReadTable(con, "v_baseline") |> 
  select(trial, ref_id, drug, arm_no, characteristic, data_type, k, n, mean, sd) |> 
  filter(characteristic %in% c(
    "Age", "Sex (n male)", "Ethnicity (n white)", "Weight", "BMI", 
    "Duration of psoriasis", "PASI", "Psoriatic arthritis", 
    "Systemic agent", "Biologic agents"
  ))
dbDisconnect(con)
