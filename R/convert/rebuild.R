## Recommend re-starting R before running to avoid errors

# Copy Access db to folder
file.copy(
  "C:/Users/p12916tc/OneDrive - The University of Manchester/Documents/@ PhD/2. Network meta-analysis/RevPal.accdb",
  "C:/Users/p12916tc/Documents/GitHub/psoriasis-rct-explorer/RevPal.accdb",
  overwrite = TRUE
)

# Convert Access db to a sqlite file
source("R/convert/convert.R")

# Generate trial estimates and CIs for forest plots
source("R/meta-analyse/trial-estimates.R")

# Load current meta-analysis results from file (does NOT re-run the 
# meta-analyses on the new data)
load("R/meta-analyse/meta-analysis.RData")

# Write meta-analysis results to sqlite file
con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

exc <- c("Izokibep", "Mirikizumab", "Phototherapy")

model_info <- bind_rows(results) |> 
  mutate(endpoint_group = if_else(likelihood == "multinomial" & 
                                    str_detect(endpoint, "pasi|dlqi"),
                                  substr(endpoint, 1, 4), endpoint)) |> 
  distinct(endpoint_group, type, likelihood, method, effects, dic) |> 
  arrange(type != "network", likelihood != "multinomial") |> 
  mutate(ma_id = row_number())

results_table <- bind_rows(results) |> 
  filter(comp_tx %notin% exc, ref_tx %notin% exc) |> 
  mutate(endpoint_group = if_else(likelihood == "multinomial" & 
                                    str_detect(endpoint, "pasi|dlqi"),
                                  substr(endpoint, 1, 4), endpoint)) |> 
  left_join(model_info, by = c("endpoint_group", "type", "likelihood", "method", "effects", "dic"),
            relationship = "many-to-one") |> 
  select(-c(type, likelihood, method, effects, dic, endpoint_group))

dbWriteTable(con, name = "ma_models", value = model_info, overwrite = TRUE)
dbWriteTable(con, name = "ma_results", value = results_table, overwrite = TRUE)

create_view_sql <- "
  CREATE VIEW v_meta_analysis AS
  SELECT *
  FROM ma_results
  LEFT JOIN ma_models ON ma_results.ma_id = ma_models.ma_id
"
dbExecute(con, "DROP VIEW IF EXISTS v_meta_analysis")
dbExecute(con, create_view_sql)

dbDisconnect(con)