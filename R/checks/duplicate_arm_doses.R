# List trials that have more than one arm with the same drug and dose.
# These arms differ by something not captured here (e.g. dosing frequency).

library(DBI)
library(RSQLite)

con <- dbConnect(SQLite(), "app/psoriasis-rcts.sqlite")

arms <- dbGetQuery(con, "
  SELECT s.trial, a.arm_no, a.arm_name, d.drug_name, a.dose_amount, du.unit_name
  FROM   arms a
  JOIN   studies s          ON s.study_id = a.study_id
  JOIN   drugs d             ON d.drug_id = a.drug_id
  LEFT   JOIN dose_units du  ON du.unit_id = a.dose_unit_id
  ORDER  BY s.trial, a.arm_no
")

dbDisconnect(con)

key <- paste(arms$trial, arms$drug_name, arms$dose_amount, arms$unit_name)
dupes <- arms[key %in% key[duplicated(key)], ]

print(dupes, row.names = FALSE)
