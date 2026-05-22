suppressPackageStartupMessages({
  library(shiny)
  library(DBI)
  library(RSQLite)
  library(DT)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

DB_PATH <- file.path(dirname(sys.frame(1)$ofile %||% "."), "revpal.sqlite")
if (!file.exists(DB_PATH)) DB_PATH <- "app/revpal.sqlite"
if (!file.exists(DB_PATH)) stop("revpal.sqlite not found - run convert.R first.")

read_db <- function(sql, params = list()) {
  con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  if (length(params)) dbGetQuery(con, sql, params = params) else dbGetQuery(con, sql)
}

drugs <- read_db(
  "SELECT DISTINCT drug FROM v_pasi WHERE drug IS NOT NULL AND drug != '' ORDER BY drug"
)$drug

# SELECT every column from one of the v_* tables, optionally filtered by a
# multi-value drug list. Each view always carries the shared id columns
# (trial, drug, dose, timepoint, timepoint_unit, n) plus its own endpoint
# columns, so callers can just `SELECT *` and trust the schema.
query_view <- function(table, drugs) {
  if (length(drugs)) {
    ph <- paste(rep("?", length(drugs)), collapse = ", ")
    read_db(
      sprintf("SELECT * FROM %s WHERE drug IN (%s)
               ORDER BY trial, arm_no, timepoint", table, ph),
      params = as.list(drugs)
    )
  } else {
    read_db(sprintf("SELECT * FROM %s ORDER BY trial, arm_no, timepoint", table))
  }
}

# "8 (10%)"; blank if responder count missing.
fmt_pasi <- function(k, n) {
  out <- rep("", length(k))
  ok  <- !is.na(k) & !is.na(n) & n > 0
  pct <- round(k[ok] / n[ok] * 100)
  out[ok] <- sprintf("%d (%d%%)", as.integer(k[ok]), pct)
  out
}

# "12.3 (4.5)"; blank if mean missing. SD optional - prints "12.3" if NA.
fmt_mean_sd <- function(mean, sd, digits = 1) {
  out <- rep("", length(mean))
  ok  <- !is.na(mean)
  m   <- formatC(mean[ok], format = "f", digits = digits)
  s   <- ifelse(is.na(sd[ok]), "",
                sprintf(" (%s)", formatC(sd[ok], format = "f", digits = digits)))
  out[ok] <- paste0(m, s)
  out
}

# "12 wks" / "3 mo" from numeric timepoint + unit code.
fmt_timepoint <- function(timepoint, unit) {
  unit_lbl <- ifelse(unit == "wk", "wks", unit)
  ifelse(is.na(timepoint), "", paste(timepoint, unit_lbl))
}

# Build a `view` reactive that pulls the right table for the current drug
# filter and runs `format_fn` (which knows how to render endpoint columns).
make_view <- function(table, format_fn) {
  function(input) {
    df <- query_view(table, input$drug)
    df$timepoint <- fmt_timepoint(df$timepoint, df$timepoint_unit)
    df$timepoint_unit <- NULL
    format_fn(df)
  }
}

# Per-view formatters: format endpoint columns and drop helper columns.
format_pasi <- function(df) {
  for (col in c("pasi50", "pasi75", "pasi90", "pasi100"))
    df[[col]] <- fmt_pasi(df[[col]], df$n)
  df[, c("trial", "drug", "dose", "timepoint", "n",
         "pasi50", "pasi75", "pasi90", "pasi100")]
}

format_pasi_abs <- function(df) {
  df$abs_pasi        <- fmt_mean_sd(df$abs_pasi_mean,        df$abs_pasi_sd)
  df$abs_pasi_change <- fmt_mean_sd(df$abs_pasi_change_mean, df$abs_pasi_change_sd)
  df[, c("trial", "drug", "dose", "timepoint", "n",
         "abs_pasi", "abs_pasi_change")]
}

format_dlqi <- function(df) {
  for (col in c("dlqi_0_1", "dlqi_0", "dlqi_5pt_dec", "dlqi_4pt_dec", "dlqi_le5"))
    df[[col]] <- fmt_pasi(df[[col]], df$n)
  df$abs_dlqi        <- fmt_mean_sd(df$abs_dlqi_mean,        df$abs_dlqi_sd)
  df$abs_dlqi_change <- fmt_mean_sd(df$abs_dlqi_change_mean, df$abs_dlqi_change_sd)
  df[, c("trial", "drug", "dose", "timepoint", "n",
         "dlqi_0_1", "dlqi_0", "dlqi_5pt_dec", "dlqi_4pt_dec", "dlqi_le5",
         "abs_dlqi", "abs_dlqi_change")]
}

format_safety <- function(df) {
  binary_cols <- c("sae", "disc_any", "disc_ae", "serious_infection",
                   "injection_site_rxn", "malignancy", "nmsc",
                   "malignancy_non_nmsc")
  for (col in binary_cols) df[[col]] <- fmt_pasi(df[[col]], df$n)
  df[, c("trial", "drug", "dose", "timepoint", "n", binary_cols)]
}

render_view <- function(df_fn, colnames) {
  renderDT({
    df <- df_fn()
    n_endpoint_cols <- ncol(df) - 5  # trial, drug, dose, timepoint, n are first 5
    datatable(
      df,
      rownames = FALSE,
      filter   = "none",
      options  = list(
        pageLength = 25,
        autoWidth  = TRUE,
        dom        = "tip",
        columnDefs = list(list(className = "dt-right",
                               targets = 3:(4 + n_endpoint_cols)))
      ),
      colnames = colnames
    )
  })
}

ui <- fluidPage(
  titlePanel("RevPal endpoint explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput("drug", "Drug",
                     choices  = drugs,
                     selected = NULL,
                     multiple = TRUE,
                     options  = list(placeholder = "(all drugs)")),
      helpText("One row per study arm × timepoint. Binary cells show",
               "n (% of arm N); continuous cells show mean (SD).")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id   = "view",
        type = "tabs",
        tabPanel("PASI thresholds", DTOutput("tbl_pasi")),
        tabPanel("Absolute PASI",   DTOutput("tbl_abs")),
        tabPanel("DLQI",            DTOutput("tbl_dlqi")),
        tabPanel("Safety",          DTOutput("tbl_safety"))
      )
    )
  )
)

server <- function(input, output, session) {
  data_pasi   <- reactive(make_view("v_pasi",     format_pasi)(input))
  data_abs    <- reactive(make_view("v_pasi_abs", format_pasi_abs)(input))
  data_dlqi   <- reactive(make_view("v_dlqi",     format_dlqi)(input))
  data_safety <- reactive(make_view("v_safety",   format_safety)(input))

  output$tbl_pasi <- render_view(data_pasi, c(
    "Trial", "Drug", "Dose", "Timepoint", "N",
    "PASI 50", "PASI 75", "PASI 90", "PASI 100"))

  output$tbl_abs <- render_view(data_abs, c(
    "Trial", "Drug", "Dose", "Timepoint", "N",
    "Absolute PASI", "Δ from baseline"))

  output$tbl_dlqi <- render_view(data_dlqi, c(
    "Trial", "Drug", "Dose", "Timepoint", "N",
    "DLQI 0/1", "DLQI 0", "5+ pt decrease", "4+ pt decrease", "DLQI ≤5",
    "Absolute DLQI", "Δ from baseline"))

  output$tbl_safety <- render_view(data_safety, c(
    "Trial", "Drug", "Dose", "Timepoint", "N",
    "Any SAE", "Disc. (any)", "Disc. (AE)", "Serious infection",
    "Injection site rxn", "Malignancy", "NMSC", "Malignancy (non-NMSC)"))
}

shinyApp(ui, server)
