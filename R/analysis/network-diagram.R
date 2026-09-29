# Ad hoc network diagram, decoupled from the Shiny app - edit REF_IDS below
# and run. Nodes = drugs (sized by total randomised patients across the
# selected trials), edges = drugs compared head-to-head within a trial
# (width = number of trials). Duplicates the app's node/edge construction
# rather than sourcing app.R, so this stays a plain standalone script.

library(DBI)
library(RSQLite)
library(visNetwork)

# Trials to include, by ref_id (= studies.study_id). Leave as NULL to include
# every trial in the database. Example below: ACCEPT 2010, AFFIRM 2022,
# AMAGINE-1 2016.
REF_IDS <- c(1, 18, 20, 24, 25, 38, 39, 42, 47, 53, 60, 71, 76, 97, 99, 105, 134, 137, 140, 144, 148, 150, 160, 162, 178, 204, 210, 217, 230, 256, 261, 272, 280, 293, 300, 303, 350, 352, 364, 382, 417, 422, 424, 439, 459, 475, 481, 489, 491, 492)

DRUG_CLASS <- c(
  "Acitretin"           = "conventional",
  "Cyclosporin"         = "conventional",
  "Fumaric acid esters" = "conventional",
  "Methotrexate"        = "conventional",
  "Phototherapy"        = "conventional",
  "Topical"             = "conventional",
  "Apremilast"          = "targeted small molecule",
  "Roflumilast"         = "targeted small molecule",
  "Deucravacitinib"     = "targeted small molecule",
  "Orismilast"          = "targeted small molecule",
  "Tofacitinib"         = "targeted small molecule",
  "Zasocitinib"         = "targeted small molecule",
  "Guselkumab"          = "il23",
  "Icotrokinra"         = "il23",
  "Mirikizumab"         = "il23",
  "Risankizumab"        = "il23",
  "Tildrakizumab"       = "il23",
  "Bimekizumab"         = "il17",
  "Brodalumab"          = "il17",
  "Ixekizumab"          = "il17",
  "Izokibep"            = "il17",
  "Netakimab"           = "il17",
  "Secukinumab"         = "il17",
  "Sonelokimab"         = "il17",
  "Xeligekimab"         = "il17",
  "Adalimumab"          = "tnf",
  "Certolizumab"        = "tnf",
  "Etanercept"          = "tnf",
  "Infliximab"          = "tnf",
  "Ustekinumab"         = "il12_23",
  "Placebo"             = "placebo"
)

CLASS_COLORS <- c(
  "il17"                    = "#E63946",
  "il23"                    = "#F4A261",
  "il12_23"                 = "#E9C46A",
  "tnf"                     = "#2A9D8F",
  "targeted small molecule" = "#264653",
  "conventional"            = "#6C757D",
  "placebo"                 = "#1D3557"
)

# Ring order for the circular layout below (Placebo sits at the centre,
# everything else arcs around it grouped by class).
CLASS_ARC_ORDER <- c("il17", "il23", "il12_23", "tnf",
                     "targeted small molecule", "conventional")

NODE_SIZE_MIN <- 14
NODE_SIZE_MAX <- 55

# --- Pull per-arm drug + patient counts for the selected trials ------------

con <- dbConnect(SQLite(), "app/psoriasis-rcts.sqlite", flags = SQLITE_RO)
ref_filter <- if (is.null(REF_IDS)) "" else
  sprintf("WHERE s.study_id IN (%s)", paste(REF_IDS, collapse = ", "))
td <- dbGetQuery(con, sprintf(
  "SELECT s.study_id AS ref_id, s.trial AS trial, a.arm_no AS arm_no,
          dr.drug_name AS drug, MAX(m.n) AS n_arm
   FROM   arms a
   JOIN   studies s        ON s.study_id = a.study_id
   LEFT   JOIN drugs dr    ON dr.drug_id = a.drug_id
   LEFT   JOIN measurements m ON m.arm_id = a.arm_id
   %s
   GROUP  BY s.study_id, s.trial, a.arm_no, dr.drug_name",
  ref_filter
))
dbDisconnect(con)

td <- td[!is.na(td$drug) & nzchar(td$drug), , drop = FALSE]

# --- Nodes: one per drug, sized by total patients ---------------------------

drugs <- sort(unique(td$drug))
n_patients <- vapply(drugs, function(d)
  sum(td$n_arm[td$drug == d], na.rm = TRUE), numeric(1))
n_trials_per_drug <- vapply(drugs, function(d)
  length(unique(td$trial[td$drug == d])), integer(1))

size_ref <- max(n_patients, 1, na.rm = TRUE)
size_from_patients <- function(n) {
  r <- sqrt(pmax(n, 1)) / sqrt(size_ref)
  pmin(NODE_SIZE_MAX, pmax(NODE_SIZE_MIN,
                           NODE_SIZE_MIN + (NODE_SIZE_MAX - NODE_SIZE_MIN) * r))
}

drug_classes <- DRUG_CLASS[drugs]
drug_classes[is.na(drug_classes)] <- "conventional"

nodes_df <- data.frame(
  id    = drugs,
  label = drugs,
  size  = size_from_patients(n_patients),
  title = sprintf("<b>%s</b><br/>%d trial(s), %s patient(s)",
                  drugs, n_trials_per_drug,
                  formatC(n_patients, format = "d", big.mark = ",")),
  color.background = unname(CLASS_COLORS[drug_classes]),
  color.border     = unname(CLASS_COLORS[drug_classes]),
  stringsAsFactors = FALSE
)

# --- Circular layout: drugs arced by class, Placebo at the centre ----------

is_placebo  <- nodes_df$id == "Placebo"
ring_drugs  <- nodes_df$id[!is_placebo]
ring_classes <- unname(drug_classes[ring_drugs])
ring_order  <- order(match(ring_classes, CLASS_ARC_ORDER), ring_drugs)
ring_drugs  <- ring_drugs[ring_order]
n_ring      <- length(ring_drugs)
angles      <- seq(0, 2 * pi, length.out = n_ring + 1)[seq_len(n_ring)]
radius      <- 500

nodes_df$x <- 0
nodes_df$y <- 0
ring_match <- match(ring_drugs, nodes_df$id)
nodes_df$x[ring_match] <- radius * cos(angles)
nodes_df$y[ring_match] <- radius * sin(angles)

# --- Edges: drug pairs sharing a trial, width = number of shared trials -----

td_pair <- unique(td[, c("trial", "drug")])
pair_rows <- list()
for (tr in unique(td_pair$trial)) {
  ds <- sort(unique(td_pair$drug[td_pair$trial == tr]))
  if (length(ds) < 2) next
  cmb <- utils::combn(ds, 2)
  pair_rows[[tr]] <- data.frame(from = cmb[1, ], to = cmb[2, ],
                                stringsAsFactors = FALSE)
}
pairs <- do.call(rbind, pair_rows)
pair_counts <- as.data.frame(table(pairs$from, pairs$to), stringsAsFactors = FALSE)
names(pair_counts) <- c("from", "to", "n_trials")
pair_counts <- pair_counts[pair_counts$n_trials > 0, ]
pair_counts$value <- pair_counts$n_trials
pair_counts$title <- sprintf("<b>%s &harr; %s</b><br/>%d trial(s)",
                             pair_counts$from, pair_counts$to,
                             pair_counts$n_trials)

# --- Render -------------------------------------------------------------

visNetwork(nodes_df, pair_counts) |>
  visNodes(shape = "dot", font = list(size = 20), physics = FALSE) |>
  visEdges(smooth  = list(enabled = FALSE),
          scaling = list(min = 2, max = 12),
          color   = list(color = "rgba(80,80,80,0.35)")) |>
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1)) |>
  visPhysics(enabled = FALSE)
