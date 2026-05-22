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

# "8 (10%)"; blank if responder count missing.
fmt_pasi <- function(k, n) {
  out <- rep("", length(k))
  ok  <- !is.na(k) & !is.na(n) & n > 0
  pct <- round(k[ok] / n[ok] * 100)
  out[ok] <- sprintf("%d (%d%%)", as.integer(k[ok]), pct)
  out
}

ui <- fluidPage(
  titlePanel("RevPal PASI explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput("drug", "Drug",
                     choices  = drugs,
                     selected = NULL,
                     multiple = TRUE,
                     options  = list(placeholder = "(all drugs)")),
      helpText("One row per study arm × timepoint. PASI cells show",
               "n (% of arm N).")
    ),
    mainPanel(
      width = 9,
      DTOutput("tbl")
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    if (length(input$drug)) {
      ph <- paste(rep("?", length(input$drug)), collapse = ", ")
      df <- read_db(
        sprintf("SELECT trial, drug, dose, timepoint, n,
                        pasi50, pasi75, pasi90, pasi100
                 FROM v_pasi
                 WHERE drug IN (%s)
                 ORDER BY trial, arm_no, timepoint", ph),
        params = as.list(input$drug)
      )
    } else {
      df <- read_db(
        "SELECT trial, drug, dose, timepoint, n,
                pasi50, pasi75, pasi90, pasi100
         FROM v_pasi
         ORDER BY trial, arm_no, timepoint"
      )
    }
    df$pasi50  <- fmt_pasi(df$pasi50,  df$n)
    df$pasi75  <- fmt_pasi(df$pasi75,  df$n)
    df$pasi90  <- fmt_pasi(df$pasi90,  df$n)
    df$pasi100 <- fmt_pasi(df$pasi100, df$n)
    df$timepoint <- ifelse(is.na(df$timepoint), "",
                           paste0(df$timepoint, " wks"))
    df
  })

  output$tbl <- renderDT({
    datatable(
      data(),
      rownames = FALSE,
      filter   = "none",
      options  = list(
        pageLength = 25,
        autoWidth  = TRUE,
        dom        = "tip",   # table + info + pagination only (no search box)
        columnDefs = list(list(className = "dt-right", targets = 3:8))
      ),
      colnames = c("Trial", "Drug", "Dose", "Timepoint",
                   "N", "PASI 50", "PASI 75", "PASI 90", "PASI 100")
    )
  })
}

shinyApp(ui, server)
