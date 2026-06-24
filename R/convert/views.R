# SQL views: the contract the Shiny app reads. Sourced by convert.R, which
# calls build_views(dst) once the relational tables are populated.
#
# One view per tab. Each carries its tab's binary responder columns AND any
# absolute / change-from-baseline columns, so the app reads a single table per
# tab:
#   v_pasi   - pasi50/75/90/100 responders, absolute & change-from-baseline
#              PASI, plus baseline PASI (outcome 11) as an arm-level column.
#   v_dlqi   - DLQI binary endpoints, absolute & change-from-baseline DLQI.
#   v_safety - binary safety outcomes.
#
# Each view is per (arm, timepoint). Columns the app expects from arm context:
# trial, ref_id, arm_no, arm_name, drug, dose, timepoint, timepoint_unit, n.

# Common arm/study context CTE. Dose is rebuilt as "<amount> <unit>"
# (integer-valued amounts render without a trailing .0).
.arm_ctx_cte <- "
WITH arm_ctx AS (
  SELECT a.arm_id,
         s.study_id     AS ref_id,
         s.trial        AS trial,
         a.arm_no       AS arm_no,
         a.arm_name     AS arm_name,
         dr.drug_name   AS drug,
         TRIM(
           CASE
             WHEN a.dose_amount IS NULL THEN ''
             WHEN a.dose_amount = CAST(a.dose_amount AS INTEGER)
               THEN CAST(CAST(a.dose_amount AS INTEGER) AS TEXT)
             ELSE CAST(a.dose_amount AS TEXT)
           END
           || ' ' || COALESCE(du.unit_name, '')
         )              AS dose,
         s.timepoint_unit AS timepoint_unit
  FROM   arms a
  JOIN   studies s          ON s.study_id = a.study_id
  LEFT   JOIN drugs dr      ON dr.drug_id = a.drug_id
  LEFT   JOIN dose_units du ON du.unit_id = a.dose_unit_id
)"

# Build one view. `pivot_cols` is a vector of pivot expressions (see case_*);
# `outcome_codes` lists the codes the view consults. `extra_join` / `extra_select`
# bolt on an arm-level LEFT JOIN (e.g. baseline PASI) that isn't part of the
# per-timepoint pivot.
.build_view <- function(dst, name, outcome_codes, pivot_cols, select_cols,
                        extra_join = "", extra_select = character(0)) {
  cat(sprintf("Building view %s ...\n", name))
  in_list <- paste(sprintf("'%s'", outcome_codes), collapse = ", ")
  pivots  <- paste(pivot_cols, collapse = ",\n           ")
  selects <- paste(c(select_cols, extra_select), collapse = ",\n           ")

  sql <- sprintf("
    CREATE VIEW %s AS
    %s,
    pivot AS (
      SELECT m.arm_id, m.timepoint,
             %s,
             MAX(m.n) AS n
      FROM   measurements m
      JOIN   outcomes o ON o.outcome_id = m.outcome_id
      WHERE  o.code IN (%s)
        AND  m.subgroup_id = 0
      GROUP BY m.arm_id, m.timepoint
    )
    SELECT ctx.trial,
           ctx.ref_id,
           ctx.arm_no,
           ctx.arm_name,
           ctx.drug,
           ctx.dose,
           pivot.timepoint AS timepoint,
           ctx.timepoint_unit,
           pivot.n         AS n,
           %s
    FROM   pivot
    JOIN   arm_ctx ctx ON ctx.arm_id = pivot.arm_id
    %s
  ", name, .arm_ctx_cte, pivots, in_list, selects, extra_join)

  dbExecute(dst, sql)
  n_rows <- dbGetQuery(dst, sprintf("SELECT COUNT(*) AS n FROM %s", name))$n
  cat(sprintf("  %s rows: %d\n", name, n_rows))
}

# MAX(CASE WHEN ... THEN <field> END) AS <alias> pivot-clause builders.
.case_k    <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.k      END) AS %s", code, alias)
.case_mean <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.mean   END) AS %s", code, alias)
.case_sd   <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.sd     END) AS %s", code, alias)
.case_med  <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.median END) AS %s", code, alias)
.case_lo   <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.lo_iqr END) AS %s", code, alias)
.case_hi   <- function(code, alias) sprintf("MAX(CASE WHEN o.code = '%s' THEN m.hi_iqr END) AS %s", code, alias)

build_views <- function(dst) {
  # v_pasi - responders + absolute/change PASI + baseline PASI (outcome 11).
  # Baseline PASI is arm-level (no timepoint), so it rides an arm-level subquery
  # join rather than the per-timepoint pivot; it repeats across the arm's rows.
  .build_view(
    dst, "v_pasi",
    outcome_codes = c("pasi50","pasi75","pasi90","pasi100",
                      "abs_pasi","abs_pasi_change"),
    pivot_cols = c(.case_k   ("pasi50",  "pasi50"),
                   .case_k   ("pasi75",  "pasi75"),
                   .case_k   ("pasi90",  "pasi90"),
                   .case_k   ("pasi100", "pasi100"),
                   .case_mean("abs_pasi",        "abs_pasi_mean"),
                   .case_sd  ("abs_pasi",        "abs_pasi_sd"),
                   .case_med ("abs_pasi",        "abs_pasi_median"),
                   .case_lo  ("abs_pasi",        "abs_pasi_lo_iqr"),
                   .case_hi  ("abs_pasi",        "abs_pasi_hi_iqr"),
                   .case_mean("abs_pasi_change", "abs_pasi_change_mean"),
                   .case_sd  ("abs_pasi_change", "abs_pasi_change_sd"),
                   .case_med ("abs_pasi_change", "abs_pasi_change_median"),
                   .case_lo  ("abs_pasi_change", "abs_pasi_change_lo_iqr"),
                   .case_hi  ("abs_pasi_change", "abs_pasi_change_hi_iqr")),
    select_cols = c("pivot.pasi50", "pivot.pasi75", "pivot.pasi90", "pivot.pasi100",
                    "pivot.abs_pasi_mean","pivot.abs_pasi_sd","pivot.abs_pasi_median",
                    "pivot.abs_pasi_lo_iqr","pivot.abs_pasi_hi_iqr",
                    "pivot.abs_pasi_change_mean","pivot.abs_pasi_change_sd",
                    "pivot.abs_pasi_change_median",
                    "pivot.abs_pasi_change_lo_iqr","pivot.abs_pasi_change_hi_iqr"),
    extra_join = "
    LEFT JOIN (
      SELECT arm_id,
             MAX(mean) AS baseline_pasi_mean,
             MAX(sd)   AS baseline_pasi_sd
      FROM   measurements
      WHERE  outcome_id = 11 AND subgroup_id = 0
      GROUP  BY arm_id
    ) bp ON bp.arm_id = pivot.arm_id",
    extra_select = c("bp.baseline_pasi_mean", "bp.baseline_pasi_sd")
  )

  # v_dlqi - DLQI binary + continuous endpoints.
  .build_view(
    dst, "v_dlqi",
    outcome_codes = c("dlqi_0_1","dlqi_0","dlqi_5pt_dec","dlqi_4pt_dec","dlqi_le5",
                      "abs_dlqi","abs_dlqi_change"),
    pivot_cols = c(.case_k   ("dlqi_0_1",        "dlqi_0_1"),
                   .case_k   ("dlqi_0",          "dlqi_0"),
                   .case_k   ("dlqi_5pt_dec",    "dlqi_5pt_dec"),
                   .case_k   ("dlqi_4pt_dec",    "dlqi_4pt_dec"),
                   .case_k   ("dlqi_le5",        "dlqi_le5"),
                   .case_mean("abs_dlqi",        "abs_dlqi_mean"),
                   .case_sd  ("abs_dlqi",        "abs_dlqi_sd"),
                   .case_med ("abs_dlqi",        "abs_dlqi_median"),
                   .case_lo  ("abs_dlqi",        "abs_dlqi_lo_iqr"),
                   .case_hi  ("abs_dlqi",        "abs_dlqi_hi_iqr"),
                   .case_mean("abs_dlqi_change", "abs_dlqi_change_mean"),
                   .case_sd  ("abs_dlqi_change", "abs_dlqi_change_sd"),
                   .case_med ("abs_dlqi_change", "abs_dlqi_change_median"),
                   .case_lo  ("abs_dlqi_change", "abs_dlqi_change_lo_iqr"),
                   .case_hi  ("abs_dlqi_change", "abs_dlqi_change_hi_iqr")),
    select_cols = c("pivot.dlqi_0_1","pivot.dlqi_0","pivot.dlqi_5pt_dec",
                    "pivot.dlqi_4pt_dec","pivot.dlqi_le5",
                    "pivot.abs_dlqi_mean","pivot.abs_dlqi_sd","pivot.abs_dlqi_median",
                    "pivot.abs_dlqi_lo_iqr","pivot.abs_dlqi_hi_iqr",
                    "pivot.abs_dlqi_change_mean","pivot.abs_dlqi_change_sd",
                    "pivot.abs_dlqi_change_median",
                    "pivot.abs_dlqi_change_lo_iqr","pivot.abs_dlqi_change_hi_iqr")
  )

  # v_safety - binary safety outcomes.
  .build_view(
    dst, "v_safety",
    outcome_codes = c("sae","disc_any","disc_ae","serious_infection",
                      "injection_site_rxn","malignancy","nmsc","malignancy_non_nmsc"),
    pivot_cols = c(.case_k("sae",                 "sae"),
                   .case_k("disc_any",            "disc_any"),
                   .case_k("disc_ae",             "disc_ae"),
                   .case_k("serious_infection",   "serious_infection"),
                   .case_k("injection_site_rxn",  "injection_site_rxn"),
                   .case_k("malignancy",          "malignancy"),
                   .case_k("nmsc",                "nmsc"),
                   .case_k("malignancy_non_nmsc", "malignancy_non_nmsc")),
    select_cols = c("pivot.sae","pivot.disc_any","pivot.disc_ae",
                    "pivot.serious_infection","pivot.injection_site_rxn",
                    "pivot.malignancy","pivot.nmsc","pivot.malignancy_non_nmsc")
  )
}
