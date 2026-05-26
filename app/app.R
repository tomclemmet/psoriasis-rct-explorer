suppressPackageStartupMessages({
  library(shiny)
  library(DBI)
  library(RSQLite)
  library(DT)
  library(visNetwork)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

DB_PATH <- file.path(dirname(sys.frame(1)$ofile %||% "."), "psoriasis-rcts.sqlite")
if (!file.exists(DB_PATH)) DB_PATH <- "app/psoriasis-rcts.sqlite"
if (!file.exists(DB_PATH)) stop("psoriasis-rcts.sqlite not found - run convert.R first.")

read_db <- function(sql, params = list()) {
  con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  if (length(params)) dbGetQuery(con, sql, params = params) else dbGetQuery(con, sql)
}

# Build network data once at startup. Nodes = drugs (sized by trial count).
# Edges = unordered drug pairs that appear together in a trial (width = number
# of trials with that head-to-head comparison). Drives the clickable NMA
# diagram which doubles as the table filter.
build_network <- function() {
  td <- read_db(
    "SELECT DISTINCT trial, drug FROM v_pasi
     WHERE drug IS NOT NULL AND drug != '' AND trial IS NOT NULL"
  )

  trial_counts <- as.data.frame(table(td$drug), stringsAsFactors = FALSE)
  names(trial_counts) <- c("id", "n_trials")
  nodes_df <- data.frame(
    id    = trial_counts$id,
    label = trial_counts$id,
    value = trial_counts$n_trials,
    title = sprintf("<b>%s</b><br/>%d trial(s)",
                    trial_counts$id, trial_counts$n_trials),
    stringsAsFactors = FALSE
  )
  nodes_df <- nodes_df[order(-nodes_df$value, nodes_df$id), ]

  # Edges: pairs within trial, dedupe with from < to, count trials.
  pair_rows <- list()
  for (tr in unique(td$trial)) {
    ds <- sort(unique(td$drug[td$trial == tr]))
    if (length(ds) < 2) next
    cmb <- utils::combn(ds, 2)
    pair_rows[[tr]] <- data.frame(from = cmb[1, ], to = cmb[2, ],
                                  stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, pair_rows)
  pair_counts <- as.data.frame(table(pairs$from, pairs$to),
                               stringsAsFactors = FALSE)
  names(pair_counts) <- c("from", "to", "n_trials")
  pair_counts <- pair_counts[pair_counts$n_trials > 0, ]
  pair_counts$id    <- sprintf("e%d", seq_len(nrow(pair_counts)))
  pair_counts$value <- pair_counts$n_trials
  pair_counts$title <- sprintf("<b>%s &harr; %s</b><br/>%d trial(s)",
                               pair_counts$from, pair_counts$to,
                               pair_counts$n_trials)

  list(nodes = nodes_df, edges = pair_counts)
}

.network <- build_network()
nodes_df <- .network$nodes
edges_df <- .network$edges

# ref_id -> DOI lookup, built once. Only ~12% of refs have a DOI in this
# dataset (URL field is empty across the board), so most trial cells stay
# plain text.
.doi_rows <- read_db(
  "SELECT study_id AS ref_id, doi AS DOI FROM publications
   WHERE is_primary = 1 AND doi IS NOT NULL AND doi != ''"
)
doi_lookup <- setNames(.doi_rows$DOI, as.character(.doi_rows$ref_id))

# Baseline PASI recorded as a "Psoriasis characteristics" outcome
# (outcome_id 11) — used as a fallback for the Absolute PASI table when an
# arm has no week-0 abs_pasi row. Curators usually record baseline PASI as
# a baseline characteristic; week 0 of the absolute-PASI longitudinal series
# is only present for the minority of arms where that timepoint was
# explicitly extracted.
.baseline_pasi <- read_db("
  SELECT a.study_id AS ref_id, a.arm_no AS arm_no,
         MAX(m.mean) AS mean, MAX(m.sd) AS sd
  FROM   measurements m
  JOIN   arms a ON a.arm_id = m.arm_id
  WHERE  m.outcome_id = 11 AND m.subgroup_id = 0
  GROUP  BY a.study_id, a.arm_no
")
.baseline_pasi_key <- paste(.baseline_pasi$ref_id, .baseline_pasi$arm_no, sep = "|")
baseline_pasi_lookup <- list(
  mean = setNames(.baseline_pasi$mean, .baseline_pasi_key),
  sd   = setNames(.baseline_pasi$sd,   .baseline_pasi_key)
)

# Wrap trial text in an <a href="https://doi.org/..."> when a DOI is known
# for that ref_id; otherwise return the trial text unchanged. Caller must
# pass `escape = FALSE` to DT for the Trial column.
fmt_trial <- function(trial, ref_id) {
  doi <- doi_lookup[as.character(ref_id)]
  has_doi <- !is.na(doi) & nzchar(doi)
  out <- trial
  out[has_doi] <- sprintf(
    '<a href="https://doi.org/%s" target="_blank" rel="noopener">%s</a>',
    doi[has_doi], htmltools::htmlEscape(trial[has_doi])
  )
  out
}

# Filter state shape:
#   NULL                                        - no filter, show all rows
#   list(kind = "node", drug = "Adalimumab")    - single-drug filter
#   list(kind = "edge", from = "A", to = "B")   - head-to-head pair: only
#                                                 rows where drug in {A,B}
#                                                 AND trial includes both
#                                                 drugs (per v_pasi).
query_view <- function(table, state) {
  base_order <- "ORDER BY trial, arm_no, timepoint"
  if (is.null(state)) {
    return(read_db(sprintf("SELECT * FROM %s %s", table, base_order)))
  }
  if (identical(state$kind, "node")) {
    return(read_db(
      sprintf("SELECT * FROM %s WHERE drug = ? %s", table, base_order),
      params = list(state$drug)
    ))
  }
  if (identical(state$kind, "edge")) {
    return(read_db(
      sprintf("SELECT * FROM %s
               WHERE drug IN (?, ?)
                 AND trial IN (
                   SELECT trial FROM v_pasi WHERE drug = ?
                   INTERSECT
                   SELECT trial FROM v_pasi WHERE drug = ?
                 )
               %s", table, base_order),
      params = list(state$from, state$to, state$from, state$to)
    ))
  }
  read_db(sprintf("SELECT * FROM %s %s", table, base_order))
}

fmt_pasi <- function(k, n) {
  out <- rep("", length(k))
  ok  <- !is.na(k) & !is.na(n) & n > 0
  pct <- round(k[ok] / n[ok] * 100)
  out[ok] <- sprintf("%d (%d%%)", as.integer(k[ok]), pct)
  out
}

fmt_mean_sd <- function(mean, sd, digits = 1) {
  out <- rep("", length(mean))
  ok  <- !is.na(mean)
  m   <- formatC(mean[ok], format = "f", digits = digits)
  s   <- ifelse(is.na(sd[ok]), "",
                sprintf(" (%s)", formatC(sd[ok], format = "f", digits = digits)))
  out[ok] <- paste0(m, s)
  out
}

fmt_timepoint <- function(timepoint, unit) {
  unit_lbl <- ifelse(unit == "wk", "wks", unit)
  ifelse(is.na(timepoint), "", paste(timepoint, unit_lbl))
}

# Build the Drug cell text: "Adalimumab 40 mg, 16 wks". Dose/timepoint
# omitted when missing.
fmt_drug <- function(drug, dose, timepoint, unit) {
  tp_txt   <- fmt_timepoint(timepoint, unit)
  has_dose <- !is.na(dose) & nzchar(dose)
  with_dose <- ifelse(has_dose, paste0(drug, " ", dose), drug)
  ifelse(nzchar(tp_txt), paste0(with_dose, ", ", tp_txt), with_dose)
}

# Per-arm baseline: value at timepoint == 0 for the same (ref_id, arm_no),
# with optional fallback to a pre-loaded `list(mean=<named>, sd=<named>)`
# keyed by "ref_id|arm_no" when no timepoint-0 row exists. Returns vectors
# aligned with df rows.
baseline_lookup <- function(df, mean_col, sd_col, fallback = NULL) {
  key <- paste(df$ref_id, df$arm_no, sep = "|")
  is_b <- !is.na(df$timepoint) & df$timepoint == 0
  bkey <- paste(df$ref_id[is_b], df$arm_no[is_b], sep = "|")
  i    <- match(key, bkey)
  mean <- df[[mean_col]][is_b][i]
  sd   <- df[[sd_col]][is_b][i]
  if (!is.null(fallback)) {
    miss <- is.na(mean)
    mean[miss] <- unname(fallback$mean[key[miss]])
    sd[miss]   <- unname(fallback$sd[key[miss]])
  }
  list(mean = mean, sd = sd)
}

# Generic binary-subset formatter: format each column as "k (pct%)", build
# the Drug cell, drop rows with no data in any of the selected endpoints,
# return Trial/Drug/N + selected columns.
format_binary_subset <- function(df, cols) {
  for (col in cols) df[[col]] <- fmt_pasi(df[[col]], df$n)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  has_any <- Reduce(`|`, lapply(cols, function(c) nzchar(df[[c]])))
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", cols)]
}

format_pasi_response <- function(df) {
  format_binary_subset(df, c("pasi50", "pasi75", "pasi90", "pasi100"))
}

format_pasi_absolute <- function(df) {
  b <- baseline_lookup(df, "abs_pasi_mean", "abs_pasi_sd",
                       fallback = baseline_pasi_lookup)
  df$baseline        <- fmt_mean_sd(b$mean, b$sd)
  df$on_tx           <- fmt_mean_sd(df$abs_pasi_mean,        df$abs_pasi_sd)
  df$abs_pasi_change <- fmt_mean_sd(df$abs_pasi_change_mean, df$abs_pasi_change_sd)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  # Drop pure baseline rows; each follow-up row now carries its baseline.
  df <- df[is.na(df$timepoint) | df$timepoint > 0, ]
  has_any <- nzchar(df$baseline) | nzchar(df$on_tx) | nzchar(df$abs_pasi_change)
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", "baseline", "on_tx", "abs_pasi_change")]
}

format_dlqi_zero <- function(df) {
  format_binary_subset(df, c("dlqi_0_1", "dlqi_0"))
}

format_dlqi_threshold <- function(df) {
  format_binary_subset(df, c("dlqi_le5"))
}

format_dlqi_change <- function(df) {
  format_binary_subset(df, c("dlqi_5pt_dec", "dlqi_4pt_dec"))
}

format_dlqi_absolute <- function(df) {
  b <- baseline_lookup(df, "abs_dlqi_mean", "abs_dlqi_sd")
  df$baseline        <- fmt_mean_sd(b$mean, b$sd)
  df$on_tx           <- fmt_mean_sd(df$abs_dlqi_mean,        df$abs_dlqi_sd)
  df$abs_dlqi_change <- fmt_mean_sd(df$abs_dlqi_change_mean, df$abs_dlqi_change_sd)
  df$drug  <- fmt_drug(df$drug, df$dose, df$timepoint, df$timepoint_unit)
  df$trial <- fmt_trial(df$trial, df$ref_id)
  df <- df[is.na(df$timepoint) | df$timepoint > 0, ]
  has_any <- nzchar(df$baseline) | nzchar(df$on_tx) | nzchar(df$abs_dlqi_change)
  df <- df[has_any, , drop = FALSE]
  df[, c("trial", "drug", "n", "baseline", "on_tx", "abs_dlqi_change")]
}

# Endpoint catalogue. Each tab has one or more endpoint groups. A group is
# defined by the source v_* table, the format fn that shapes it for display,
# and the column headers shown to the user.
endpoint_groups <- list(
  pasi = list(
    label = "PASI",
    groups = list(
      response = list(
        label    = "PASI 50 / 75 / 90 / 100",
        table    = "v_pasi",
        fmt      = format_pasi_response,
        colnames = c("Trial", "Drug", "N",
                     "PASI 50", "PASI 75", "PASI 90", "PASI 100")
      ),
      absolute = list(
        label    = "Absolute PASI",
        table    = "v_pasi_abs",
        fmt      = format_pasi_absolute,
        colnames = c("Trial", "Drug", "N",
                     "Baseline", "Follow-up", "Δ from baseline")
      )
    )
  ),
  dlqi = list(
    label = "DLQI",
    groups = list(
      zero = list(
        label    = "DLQI 0/1, DLQI 0",
        table    = "v_dlqi",
        fmt      = format_dlqi_zero,
        colnames = c("Trial", "Drug", "N", "DLQI 0/1", "DLQI 0")
      ),
      threshold = list(
        label    = "DLQI ≤ 5",
        table    = "v_dlqi",
        fmt      = format_dlqi_threshold,
        colnames = c("Trial", "Drug", "N", "DLQI ≤ 5")
      ),
      change = list(
        label    = "5+ / 4+ point decrease",
        table    = "v_dlqi",
        fmt      = format_dlqi_change,
        colnames = c("Trial", "Drug", "N",
                     "5+ pt decrease", "4+ pt decrease")
      ),
      absolute = list(
        label    = "Absolute DLQI",
        table    = "v_dlqi",
        fmt      = format_dlqi_absolute,
        colnames = c("Trial", "Drug", "N",
                     "Baseline", "Follow-up", "Δ from baseline")
      )
    )
  ),
  safety = list(
    label = "Safety",
    groups = list(
      sae = list(
        label    = "Any SAE",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df, c("sae")),
        colnames = c("Trial", "Drug", "N", "Any SAE")
      ),
      disc = list(
        label    = "Discontinuation (any, due to AE)",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df,
                     c("disc_any", "disc_ae")),
        colnames = c("Trial", "Drug", "N",
                     "Disc. (any)", "Disc. (AE)")
      ),
      serious_infection = list(
        label    = "Serious infections",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df, c("serious_infection")),
        colnames = c("Trial", "Drug", "N", "Serious infection")
      ),
      injection_site_rxn = list(
        label    = "Injection-site reactions",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df, c("injection_site_rxn")),
        colnames = c("Trial", "Drug", "N", "Injection site rxn")
      ),
      malignancy = list(
        label    = "Malignancy, NMSC, malignancy (non-NMSC)",
        table    = "v_safety",
        fmt      = function(df) format_binary_subset(df,
                     c("malignancy", "nmsc", "malignancy_non_nmsc")),
        colnames = c("Trial", "Drug", "N",
                     "Malignancy", "NMSC", "Malignancy (non-NMSC)")
      )
    )
  )
)

ui <- fluidPage(
  tags$head(tags$script(HTML("
    // Disable specific <option> values inside a Shiny select input. Used to
    // grey out endpoint groups that have zero rows under the current filter.
    Shiny.addCustomMessageHandler('set_disabled_options', function(msg) {
      var sel = document.getElementById(msg.input_id);
      if (!sel) return;
      var disabled = Array.isArray(msg.disabled) ? msg.disabled : [];
      Array.from(sel.options).forEach(function(opt) {
        opt.disabled = disabled.indexOf(opt.value) !== -1;
      });
    });
  "))),
  tags$head(tags$style(HTML("
    .filter-bar { padding: 8px 12px; background: #f4f6f9; border-radius: 6px;
                  margin: 6px 0 12px 0; display: flex; align-items: center;
                  gap: 12px; font-size: 15px; }
    .filter-bar .label { color: #555; }
    .filter-bar .value { font-weight: 600; color: #1F4E8C; }
    .filter-bar .btn { padding: 2px 10px; font-size: 13px; }
    #nma { background: #fafbfc; border: 1px solid #e3e6ea; border-radius: 6px; }
    .endpoint-picker { margin: 0; padding: 6px 0; background: #ffffff;
                       height: 46px; box-sizing: border-box; }
    .endpoint-picker .form-group { margin-bottom: 0; }
    .endpoint-picker select { width: 100%; height: 34px; padding: 4px 8px;
                              border: 1px solid #ccc; border-radius: 4px;
                              background: #ffffff; font-size: 14px; }
    .endpoint-picker select option:disabled { color: #b0b0b0; }
    /* Use flexbox so left/right columns share the row's full height; sticky
       then has room to stick as the user scrolls the table. */
    /* Lock the page so only the table column scrolls; the diagram never moves. */
    html, body { height: 100%; overflow: hidden; }
    .container-fluid, .container { height: 100%; display: flex; flex-direction: column; }
    .row.split { display: flex; align-items: stretch; flex: 1 1 auto;
                 min-height: 0; }
    .row.split > [class*='col-'] { float: none; }
    /* Only the table body scrolls; tabs stay pinned at the top of the column. */
    .row.split .col-table { display: flex; flex-direction: column;
                            max-height: 100%; min-height: 0; }
    .row.split .col-table .tabbable { display: flex; flex-direction: column;
                                       flex: 1 1 auto; min-height: 0; }
    .row.split .col-table .tab-content { flex: 1 1 auto; overflow: auto;
                                          min-height: 0; }
    /* Pin the endpoint dropdown at the top of the scrolling tab-content,
       and pin the table header directly underneath it. */
    .row.split .col-table .tab-content .endpoint-picker {
      position: sticky; top: 0; z-index: 3;
    }
    .row.split .col-table .tab-content table.dataTable thead th {
      position: sticky; top: 46px; background: #ffffff; z-index: 2;
      box-shadow: inset 0 -1px 0 #ddd;
    }
    .app-footer { margin-top: 12px; font-size: 12px; color: #888; }
    .app-footer a { color: #1F4E8C; text-decoration: none; }
    .app-footer a:hover { text-decoration: underline; }
  "))),
  titlePanel("Psoriasis RCT Explorer"),
  div(class = "filter-bar",
      span(class = "label", "Filter:"),
      uiOutput("filter_label", inline = TRUE),
      actionButton("clear_filter", "Clear", class = "btn btn-default"),
      downloadButton("download_db", "Download SQLite",
                     class = "btn btn-default")
  ),
  fluidRow(class = "split",
    column(6,
      visNetworkOutput("nma", height = "640px"),
      helpText("Click a node to filter to one drug; click an edge to show only",
               "trials comparing that pair. Click empty space, or the Clear",
               "button, to reset."),
      tags$footer(class = "app-footer",
        HTML("&copy; 2026 Thomas Clemmet"),
        " · ",
        tags$a(href = "https://github.com/tomclemmet/psoriasis-rct-explorer",
               target = "_blank", rel = "noopener", "GitHub"),
        " · ",
        tags$a(href = "https://www.crd.york.ac.uk/PROSPERO/view/CRD420261306630",
               target = "_blank", rel = "noopener", "PROSPERO record")
      )
    ),
    column(6, class = "col-table",
      do.call(tabsetPanel, c(
        list(id = "view", type = "tabs"),
        lapply(names(endpoint_groups), function(tab_id) {
          tab <- endpoint_groups[[tab_id]]
          group_choices <- setNames(names(tab$groups),
                                    vapply(tab$groups, `[[`, "", "label"))
          tabPanel(
            tab$label,
            div(class = "endpoint-picker",
                selectInput(paste0("group_", tab_id),
                            label = NULL,
                            choices  = group_choices,
                            selected = group_choices[[1]],
                            selectize = FALSE,
                            width    = "100%")),
            DTOutput(paste0("tbl_", tab_id))
          )
        })
      ))
    )
  )
)

server <- function(input, output, session) {

  filter_state <- reactiveVal(NULL)

  output$filter_label <- renderUI({
    s <- filter_state()
    if (is.null(s))                  span(class = "value", "all drugs")
    else if (s$kind == "node")       span(class = "value", s$drug)
    else if (s$kind == "edge")       span(class = "value",
                                          sprintf("%s ↔ %s (head-to-head)",
                                                  s$from, s$to))
  })

  # For each (tab, group), does the current filter yield any rows after
  # formatting? Cache table queries within one pass since v_safety is reused
  # across five groups.
  availability <- reactive({
    state <- filter_state()
    cache <- list()
    get_tbl <- function(tbl) {
      if (is.null(cache[[tbl]])) cache[[tbl]] <<- query_view(tbl, state)
      cache[[tbl]]
    }
    lapply(endpoint_groups, function(tab) {
      vapply(tab$groups, function(grp) {
        n <- tryCatch(nrow(grp$fmt(get_tbl(grp$table))),
                      error = function(e) 0L)
        isTRUE(n > 0)
      }, logical(1))
    })
  })

  # Push disabled-option lists into each tab's native <select>. We don't
  # auto-switch the selection: that would re-render the same DT widget with
  # a different column count and DataTables.js fires a "column not found"
  # warning. Keeping the user's pick lets the empty-table path (which has
  # always worked) handle the no-data case.
  observe({
    av <- availability()
    for (tab_id in names(av)) {
      avail <- av[[tab_id]]
      session$sendCustomMessage("set_disabled_options", list(
        input_id = paste0("group_", tab_id),
        # Wrap with I() so an empty vector serializes as `[]` (a bare empty
        # list becomes `{}`, which would break indexOf in the JS handler).
        disabled = I(as.character(names(avail)[!avail]))
      ))
    }
  })

  # Build one renderDT per tab. Reactive picks the active endpoint group
  # from that tab's dropdown, queries the right table, formats it.
  for (tab_id in names(endpoint_groups)) local({
    this_tab    <- tab_id
    tab_cfg     <- endpoint_groups[[this_tab]]
    output[[paste0("tbl_", this_tab)]] <- renderDT({
      gid <- input[[paste0("group_", this_tab)]]
      req(gid)
      grp <- tab_cfg$groups[[gid]]
      df <- grp$fmt(query_view(grp$table, filter_state()))
      n_endpoint_cols <- ncol(df) - 3
      datatable(
        df,
        rownames = FALSE,
        filter   = "none",
        escape   = -1,  # Only the Trial column (col 1) contains HTML (an
                        # <a href="https://doi.org/..."> link when a DOI is
                        # known); fmt_trial() escapes the visible trial name
                        # itself. Every other column stays escaped.
        options  = list(
          pageLength = 25,
          autoWidth  = FALSE,
          dom        = "tip",
          columnDefs = list(list(className = "dt-right",
                                 targets = 2:(2 + n_endpoint_cols))),
          # On page change, scroll the surrounding .tab-content (our custom
          # scroll container) back to the top — DT only scrolls its own
          # internal viewport, which we don't use.
          initComplete = JS(
            "function() {",
            "  var api = this.api();",
            "  api.on('page.dt', function() {",
            "    var el = $(api.table().node()).closest('.tab-content')[0];",
            "    if (el) el.scrollTo({ top: 0 });",
            "  });",
            "}"
          )
        ),
        colnames = grp$colnames
      )
    })
  })

  output$nma <- renderVisNetwork({
    visNetwork(nodes_df, edges_df) |>
      visIgraphLayout(layout    = "layout_with_kk",
                      randomSeed = 42,
                      physics    = FALSE) |>
      visNodes(shape   = "dot",
               scaling = list(min = 18, max = 55,
                              label = list(enabled = TRUE,
                                           min = 22, max = 34)),
               font    = list(size = 28, face = "Helvetica",
                              strokeWidth = 4, strokeColor = "#ffffff"),
               color   = list(background = "#4C9AFF",
                              border     = "#1F4E8C",
                              highlight  = list(background = "#FF8A3D",
                                                border     = "#B5521A"),
                              hover      = list(background = "#7FB5FF",
                                                border     = "#1F4E8C")),
               borderWidth = 2) |>
      visEdges(smooth   = list(enabled = TRUE, type = "continuous"),
               scaling  = list(min = 1, max = 10),
               color    = list(color = "rgba(80,80,80,0.30)",
                               highlight = "#FF8A3D",
                               hover     = "#1F4E8C")) |>
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1,
                                         hover = TRUE,
                                         labelOnly = FALSE),
                 nodesIdSelection = FALSE) |>
      visInteraction(navigationButtons = FALSE, multiselect = FALSE,
                     tooltipDelay = 150, hover = TRUE,
                     zoomView = TRUE, dragView = TRUE) |>
      visEvents(
        selectNode   = "function(p){ Shiny.setInputValue('nma_node', p.nodes[0], {priority:'event'}); }",
        selectEdge   = "function(p){ if(p.nodes && p.nodes.length) return;
                                     Shiny.setInputValue('nma_edge', p.edges[0], {priority:'event'}); }",
        deselectNode = "function(p){ Shiny.setInputValue('nma_clear', Math.random(), {priority:'event'}); }",
        deselectEdge = "function(p){ Shiny.setInputValue('nma_clear', Math.random(), {priority:'event'}); }"
      )
  })

  # Node click -> single-drug filter
  observeEvent(input$nma_node, {
    req(input$nma_node)
    filter_state(list(kind = "node", drug = input$nma_node))
  })

  # Edge click -> head-to-head pair filter
  observeEvent(input$nma_edge, {
    req(input$nma_edge)
    e <- edges_df[edges_df$id == input$nma_edge, , drop = FALSE]
    if (!nrow(e)) return()
    filter_state(list(kind = "edge", from = e$from[1], to = e$to[1]))
  })

  # Empty-canvas click -> clear
  observeEvent(input$nma_clear, { filter_state(NULL) })

  output$download_db <- downloadHandler(
    filename    = function() "psoriasis-rcts.sqlite",
    contentType = "application/x-sqlite3",
    content     = function(file) file.copy(DB_PATH, file)
  )
  observeEvent(input$clear_filter, {
    filter_state(NULL)
    visNetworkProxy("nma") |> visUnselectAll()
  })
}

shinyApp(ui, server)
