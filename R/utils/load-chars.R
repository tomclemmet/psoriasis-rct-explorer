# Loads relevant baseline characteristics
library(DBI)
library(dplyr)

# Load data
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

my_chars <- c(
  "Age", "Sex (n male)", "Weight",  "PASI", "Biologic agents"
)

chars <- dbReadTable(con, "v_baseline") |> 
  select(trial, ref_id, drug, arm_no, characteristic, data_type, k, n, mean, sd) |> 
  filter(characteristic %in% my_chars) |>
  mutate(.by = ref_id, arm_no = dense_rank(arm_no)) |> 
  filter(ref_id %notin% reg_ids) |> 
  filter(ref_id %notin% c(16, 128, 132, 167, 170, 181, 188, 239, 283, 287, 315, 369, 460))

dbDisconnect(con)

reg_ids <- chars |> 
  filter(ref_id %notin% c(16, 128, 132, 167, 170, 181, 188, 239, 283, 287, 315, 369, 460)) |> 
  distinct(ref_id, characteristic) |> 
  group_by(ref_id) |> 
  filter(any(characteristic == "Biologic agents")) |> 
  filter(any(characteristic == "Weight")) |> 
  filter(any(characteristic == "Age")) |> 
  filter(any(characteristic == "PASI")) |> 
  # filter(any(characteristic == "Ethnicity (n white)")) |>
  # filter(any(characteristic == "Duration of psoriasis")) |> 
  filter(any(characteristic == "Sex (n male)")) |> 
  distinct(ref_id) |> pull(ref_id)

char_list <- list()

for (i in 1:length(my_chars)) {
  
  char_list[[my_chars[i]]] <- chars |> 
    mutate(x = if_else(data_type == "Continuous", mean, k/n)) |> 
    filter(characteristic == my_chars[i], !is.na(x)) |> 
    select(ref_id, arm_no, x) |> 
    pivot_wider(names_from = arm_no, values_from = x, names_sort = TRUE) |> 
    select(-ref_id) |> 
    as.matrix()
  
}
names(char_list) <- c("age", "sex", "weight", "pasi", "bio")
