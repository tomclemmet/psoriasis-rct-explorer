fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f (%.*f to %.*f)",
          digits, est, digits, lo, digits, hi)
}

dic_extra <- function(dic) if (length(dic) && !is.na(dic)) c(DIC = sprintf("%.1f", dic)) else NULL

ma_tooltip <- function(label, est, lo, hi, extra = NULL, digits = 2) {
  parts <- c(
    label,
    sprintf("Effect: %s", fmt_ci(est, lo, hi, digits))
  )
  if (length(extra))
    parts <- c(parts, sprintf("%s: %s", names(extra), as.character(extra)))
  paste(parts, collapse = "\n")
}

forest_ticks <- function(scale, xmin, xmax) {
  if (identical(scale, "rr")) {
    candidates <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
    candidates[candidates >= xmin & candidates <= xmax]
  } else if (identical(scale, "prop")) {
    span <- xmax - xmin
    step <- if      (span <= 0.10) 0.02
            else if (span <= 0.20) 0.05
            else if (span <= 0.50) 0.10
            else                   0.25
    seq(0, xmax, by = step)
  } else {
    pretty(c(xmin, xmax), n = 5)
  }
}

forest_tick_label <- function(scale, v) {
  if (identical(scale, "rr"))    sub("\\.0$", "", formatC(v, format = "g"))
  else if (identical(scale, "prop")) sprintf("%g", v)
  else                            formatC(v, format = "g")
}

forest_xscale <- function(scale, xmin, xmax, plot_left, plot_w) {
  if (identical(scale, "rr")) {
    lmin <- log10(xmin); lmax <- log10(xmax)
    function(x) plot_left + (log10(pmax(x, xmin / 10)) - lmin) /
                            (lmax - lmin) * plot_w
  } else {
    function(x) plot_left + (x - xmin) / (xmax - xmin) * plot_w
  }
}

forest_xlimits <- function(scale, rows, pooled) {
  vals <- c(rows$lo, rows$hi)
  if (!is.null(pooled) && nrow(pooled)) vals <- c(vals, pooled$lo, pooled$hi)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(switch(scale, rr = c(0.1, 10), prop = c(0, 1),
                                   c(-1, 1)))
  if (identical(scale, "rr")) {
    candidates <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
    lo <- min(vals, 1); hi <- max(vals, 1)
    lr  <- log10(hi) - log10(lo)
    pad <- max(lr * 0.05, 0.05)
    lo  <- 10 ^ (log10(lo) - pad)
    hi  <- 10 ^ (log10(hi) + pad)
    lo <- max(candidates[candidates <= lo], 0.01)
    hi <- min(candidates[candidates >= hi], 200)
    c(lo, hi)
  } else if (identical(scale, "prop")) {
    hi_data <- max(vals)
    if (hi_data <= 0.5) {
      pad <- max(hi_data * 0.10, 0.02)
      hi_snap <- ceiling((hi_data + pad) / 0.05) * 0.05
      c(0, min(hi_snap, 0.5))
    } else {
      c(0, 1)
    }
  } else {
    pad <- (max(vals) - min(vals)) * 0.1
    c(min(min(vals), 0) - pad, max(max(vals), 0) + pad)
  }
}

forest_diamond_points <- function(xL, xC, xR, yC, hh) {
  sprintf("%g,%g %g,%g %g,%g %g,%g",
          xL, yC, xC, yC - hh, xR, yC, xC, yC + hh)
}

# Pooled diamonds are grouped under one heading per (type, <first
# MA_AXIS_COLUMNS entry>) — e.g. "Multinomial network estimates" — with each
# diamond below labelled by only the *remaining* axes + effects (e.g.
# "Standard, FE" / "Class effects, RE"). This keeps on-chart text short even
# as more axis values accumulate; the full detail still lives in the
# tooltip. MA_AXIS_COLUMNS is defined in app.R (visible here since forest.R
# is source()'d into app.R's environment) — MA_AXIS_COLUMNS[1] (likelihood
# today) is the heading axis, the rest (method today) are row axes.
cap1 <- function(s) { substr(s, 1, 1) <- toupper(substr(s, 1, 1)); s }

# Heading text for a block of pooled diamonds sharing (type, heading_value).
# `heading_value` is NA when that axis doesn't vary among this endpoint's
# combinations, so the heading stays plain ("Network estimates") until
# there's actually more than one value to disambiguate ("Multinomial network
# estimates").
pooled_group_heading <- function(type, heading_value) {
  scope_lc <- switch(type, network = "network", pairwise = "pairwise",
                     univariate = "pooled", tolower(as.character(type)))
  prefix <- if (!is.null(heading_value) && !is.na(heading_value))
    paste0(tolower(axis_display_label(heading_value)), " ") else ""
  cap1(sprintf("%s%s estimates", prefix, scope_lc))
}

# CSS class + short row label + full tooltip header for one pooled diamond
# within a pooled_group_heading() block. `row_axis_values` is a named list of
# every MA_AXIS_COLUMNS entry except the heading axis (method today); an
# entry is NA when that axis doesn't vary for this block, so the row label
# stays as plain as possible ("FE" alone, or "Class effects, FE" once method
# is also ambiguous). The tooltip header always spells out every axis in
# full, regardless of what's collapsed into the heading on-chart.
pooled_row_meta <- function(type, effects, heading_value, row_axis_values) {
  eff_abbrev <- switch(effects, fixed = "FE", random = "RE",
                       toupper(substr(as.character(effects), 1, 2)))
  eff_word   <- switch(effects, fixed = "Fixed effects", random = "Random effects",
                       tools::toTitleCase(as.character(effects)))

  row_words <- unlist(lapply(row_axis_values, function(v) {
    if (is.null(v) || !length(v) || is.na(v)) NULL else cap1(tolower(axis_display_label(v)))
  }), use.names = FALSE)
  label <- paste(c(row_words, eff_abbrev), collapse = ", ")

  scope_lc <- switch(type, network = "network", pairwise = "pairwise",
                     univariate = "pooled", tolower(as.character(type)))
  axis_words <- unlist(lapply(c(list(heading_value), row_axis_values), function(v) {
    if (is.null(v) || !length(v) || is.na(v)) NULL else tolower(axis_display_label(v))
  }), use.names = FALSE)
  phrase <- paste(c(axis_words, scope_lc, "estimate"), collapse = " ")

  list(
    klass  = if (identical(effects, "random")) "ma-pooled-re" else "ma-pooled-fe",
    label  = label,
    tt_hdr = sprintf("%s %s", eff_word, phrase)
  )
}

forest_svg <- function(rows, pooled, scale = "rr", width = 880,
                       axis_label = NULL, dir_left = NULL, dir_right = NULL) {
  esc <- function(x) htmltools::htmlEscape(x, attribute = FALSE)
  esc_attr <- function(x) htmltools::htmlEscape(x, attribute = TRUE)
  if (!nrow(rows) && (is.null(pooled) || !nrow(pooled))) {
    return(sprintf('<div class="ma-empty">No meta-analysable data for this comparison.</div>'))
  }

  ROW_H            <- 22
  HEADER_PAD       <- 8
  POOLED_GAP       <- 14
  POOLED_H         <- 22
  POOLED_HEADING_H <- 20
  AXIS_PAD         <- 28
  AXIS_LABEL_H <- if (length(axis_label) || length(dir_left) || length(dir_right)) 50 else 0
  BOTTOM_PAD   <- 6
  LEFT_MARGIN  <- 220
  RIGHT_MARGIN <- 175
  PLOT_LEFT    <- LEFT_MARGIN
  PLOT_W       <- width - LEFT_MARGIN - RIGHT_MARGIN

  n_trials  <- nrow(rows)
  n_pooled  <- if (!is.null(pooled)) nrow(pooled) else 0L

  # Sort pooled rows into heading-groups: by type (pairwise/univariate
  # sections before network), then the heading axis (MA_AXIS_COLUMNS[1],
  # likelihood today) alphabetically, then the remaining row axes (default
  # value first, e.g. "standard" leads "class_effects"), then effects (fixed
  # before random) — and tag which row starts a new heading. The heading axis
  # sorts plain-alphabetically rather than default-first so network sections
  # read "Binomial, Multinomial" — a fixed display convention, unlike the
  # row/toggle axes which prefer MA_AXIS_DEFAULT.
  ma_heading_col <- MA_AXIS_COLUMNS[1]
  ma_row_cols    <- if (length(MA_AXIS_COLUMNS) > 1) MA_AXIS_COLUMNS[-1] else character(0)
  n_pooled_groups <- 0L
  pooled_group_start <- logical(0)
  if (n_pooled) {
    axis_sort_rank <- function(values, ordered_vals) {
      ifelse(is.na(values), 0L, match(values, ordered_vals))
    }
    type_rank <- c(pairwise = 1L, univariate = 1L, network = 2L)
    heading_ordered <- sort(unique(pooled[[ma_heading_col]][!is.na(pooled[[ma_heading_col]])]))
    sort_keys <- list(
      unname(type_rank[pooled$type]),
      axis_sort_rank(pooled[[ma_heading_col]], heading_ordered)
    )
    for (rc in ma_row_cols) {
      row_ordered <- order_axis_values(rc, pooled[[rc]][!is.na(pooled[[rc]])])
      sort_keys[[length(sort_keys) + 1L]] <- axis_sort_rank(pooled[[rc]], row_ordered)
    }
    sort_keys[[length(sort_keys) + 1L]] <- match(pooled$effects, c("fixed", "random"))
    pooled <- pooled[do.call(order, sort_keys), , drop = FALSE]

    group_key <- paste(pooled$type, ifelse(is.na(pooled[[ma_heading_col]]), "", pooled[[ma_heading_col]]))
    pooled_group_start <- c(TRUE, group_key[-1] != group_key[-length(group_key)])
    n_pooled_groups <- sum(pooled_group_start)
  }

  row_heights <- if (!is.null(rows$row_h))     rows$row_h     else rep(ROW_H, n_trials)
  gaps_above  <- if (!is.null(rows$gap_above)) rows$gap_above else rep(0L,    n_trials)

  trials_body_h <- sum(row_heights) + sum(gaps_above)
  body_h    <- trials_body_h +
    (if (n_pooled) POOLED_GAP + n_pooled_groups * POOLED_HEADING_H + n_pooled * POOLED_H else 0)
  height    <- HEADER_PAD + body_h + AXIS_PAD + AXIS_LABEL_H + BOTTOM_PAD

  xlim  <- forest_xlimits(scale, rows, pooled)
  xmin  <- xlim[1]; xmax <- xlim[2]
  xfn   <- forest_xscale(scale, xmin, xmax, PLOT_LEFT, PLOT_W)
  ticks <- forest_ticks(scale, xmin, xmax)

  ref_line <- switch(scale, rr = 1, md = 0, prop = NA_real_)
  plot_bottom <- HEADER_PAD + body_h
  plot_top    <- HEADER_PAD - 4

  yc_vals <- numeric(n_trials)
  y_run <- HEADER_PAD
  for (.i in seq_len(n_trials)) {
    y_run       <- y_run + gaps_above[.i]
    yc_vals[.i] <- y_run + row_heights[.i] / 2
    y_run       <- y_run + row_heights[.i]
  }

  ci_digits <- if (identical(scale, "prop")) 3 else 2

  square_n <- if (!is.null(rows$square_n)) rows$square_n else rep(NA_real_, n_trials)
  max_n <- suppressWarnings(max(square_n, na.rm = TRUE))
  square_size <- if (!is.finite(max_n) || max_n <= 0) {
    rep(8, n_trials)
  } else {
    pmax(5, pmin(12, 5 + 7 * sqrt(square_n / max_n)))
  }
  square_size[is.na(square_size)] <- 8

  parts <- character(0)

  parts[length(parts) + 1L] <- sprintf(
    '<svg class="ma-forest" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%" height="%dpx" preserveAspectRatio="xMinYMin meet" role="img">',
    width, as.integer(height), as.integer(height))

  parts[length(parts) + 1L] <- sprintf(
    '<rect class="ma-plot-bg" x="%g" y="%g" width="%g" height="%g"/>',
    PLOT_LEFT, plot_top, PLOT_W, plot_bottom - plot_top)

  if (!is.na(ref_line) && ref_line >= xmin && ref_line <= xmax) {
    rx <- xfn(ref_line)
    parts[length(parts) + 1L] <- sprintf(
      '<line class="ma-refline" x1="%g" y1="%g" x2="%g" y2="%g"/>',
      rx, plot_top, rx, plot_bottom)
  }

  X_MIN_PX <- PLOT_LEFT
  X_MAX_PX <- PLOT_LEFT + PLOT_W

  emit_data_row <- function(yc, est, lo, hi, sz, label_left, ci_text,
                            tooltip, klass = "ma-square-default",
                            row_h = ROW_H, badge = "") {
    xL_raw <- xfn(lo); xR_raw <- xfn(hi); xE_raw <- xfn(est)
    finite_or <- function(x, fallback) if (is.finite(x)) x else fallback
    xL <- max(X_MIN_PX, finite_or(xL_raw, X_MIN_PX))
    xR <- min(X_MAX_PX, finite_or(xR_raw, X_MAX_PX))
    xE_raw_safe <- finite_or(xE_raw, (X_MIN_PX + X_MAX_PX) / 2)
    xE_clamped <- min(max(xE_raw_safe, X_MIN_PX + sz / 2), X_MAX_PX - sz / 2)
    off_left  <- is.finite(xL_raw) && xL_raw < X_MIN_PX
    off_right <- is.finite(xR_raw) && xR_raw > X_MAX_PX
    has_data  <- is.finite(xE_raw)
    label_x <- if (nzchar(badge)) LEFT_MARGIN - 30 else LEFT_MARGIN - 10
    paste0(
      '<g class="ma-row" data-tt="', esc_attr(tooltip), '">',
      sprintf('<text class="ma-rowlabel" x="%g" y="%g">%s</text>',
              label_x, yc + 4, esc(label_left)),
      if (nzchar(badge))
        sprintf('<text class="ma-rowlabel-sub" x="%g" y="%g">%s</text>',
                LEFT_MARGIN - 10, yc + 4, esc(badge))
      else "",
      sprintf('<rect class="ma-rowhit" x="%g" y="%g" width="%g" height="%g"/>',
              PLOT_LEFT, yc - row_h / 2 + 1, PLOT_W, row_h - 2),
      if (has_data && xL < xR)
        sprintf('<line class="ma-ci" x1="%g" y1="%g" x2="%g" y2="%g"/>',
                xL, yc, xR, yc) else "",
      if (off_left)
        sprintf('<polygon class="ma-ci-arrow" points="%g,%g %g,%g %g,%g"/>',
                X_MIN_PX, yc, X_MIN_PX + 6, yc - 4, X_MIN_PX + 6, yc + 4)
      else "",
      if (off_right)
        sprintf('<polygon class="ma-ci-arrow" points="%g,%g %g,%g %g,%g"/>',
                X_MAX_PX, yc, X_MAX_PX - 6, yc - 4, X_MAX_PX - 6, yc + 4)
      else "",
      if (has_data)
        sprintf('<rect class="ma-square %s" x="%g" y="%g" width="%g" height="%g"/>',
                klass, xE_clamped - sz / 2, yc - sz / 2, sz, sz)
      else "",
      sprintf('<text class="ma-citext" x="%g" y="%g">%s</text>',
              width - RIGHT_MARGIN + 10, yc + 4, esc(ci_text)),
      '</g>'
    )
  }

  row_klass <- if (!is.null(rows$klass)) rows$klass
               else rep("ma-square-default", n_trials)

  row_badge <- if (!is.null(rows$badge)) rows$badge else rep("", n_trials)
  for (i in seq_len(n_trials)) {
    parts[length(parts) + 1L] <- emit_data_row(
      yc_vals[i], rows$est[i], rows$lo[i], rows$hi[i],
      sz          = square_size[i],
      label_left  = rows$label[i],
      ci_text     = fmt_ci(rows$est[i], rows$lo[i], rows$hi[i], ci_digits),
      tooltip     = rows$tooltip[i] %||% rows$label[i],
      klass       = row_klass[i],
      row_h       = row_heights[i],
      badge       = row_badge[i]
    )
  }

  if (n_pooled) {
    y_run <- HEADER_PAD + trials_body_h + POOLED_GAP
    for (i in seq_len(n_pooled)) {
      # Styling/label derived from the (type, effects, MA_AXIS_COLUMNS...)
      # columns rather than a hardcoded per-combination lookup, so a pooled
      # row for a new axis value (or a new `type`) renders sensibly with zero
      # code changes here — see pooled_group_heading()/pooled_row_meta() above
      # and pooled_diamonds() in app.R, which is the sole place `pooled` rows
      # get built.
      type_i    <- if (!is.null(pooled$type))    pooled$type[i]    else NA_character_
      effects_i <- if (!is.null(pooled$effects)) pooled$effects[i] else NA_character_
      heading_i <- if (!is.null(pooled[[ma_heading_col]])) pooled[[ma_heading_col]][i] else NA_character_

      if (pooled_group_start[i]) {
        y_run <- y_run + POOLED_HEADING_H
        parts[length(parts) + 1L] <- sprintf(
          '<text class="ma-pooled-heading" x="%g" y="%g">%s</text>',
          LEFT_MARGIN - 10, y_run - POOLED_HEADING_H / 2 + 4,
          esc(pooled_group_heading(type_i, heading_i)))
      }
      yc <- y_run + POOLED_H / 2

      xL_raw <- xfn(pooled$lo[i]); xR_raw <- xfn(pooled$hi[i])
      xE_raw <- xfn(pooled$est[i])
      xL <- max(X_MIN_PX, if (is.finite(xL_raw)) xL_raw else X_MIN_PX)
      xR <- min(X_MAX_PX, if (is.finite(xR_raw)) xR_raw else X_MAX_PX)
      xE <- if (is.finite(xE_raw)) min(max(xE_raw, X_MIN_PX), X_MAX_PX)
            else (X_MIN_PX + X_MAX_PX) / 2

      row_axis_i <- setNames(lapply(ma_row_cols, function(col)
        if (!is.null(pooled[[col]])) pooled[[col]][i] else NA_character_), ma_row_cols)
      meta   <- pooled_row_meta(type_i, effects_i, heading_i, row_axis_i)
      klass    <- meta$klass
      kind_lbl <- meta$label
      pts    <- forest_diamond_points(xL, xE, xR, yc, 7)
      ci_str <- fmt_ci(pooled$est[i], pooled$lo[i], pooled$hi[i], ci_digits)
      tt_hdr <- meta$tt_hdr
      dic_i <- if (!is.null(pooled$dic)) pooled$dic[i] else NA_real_
      tt <- if (!is.na(dic_i))
        sprintf("%s\n%s\nDIC: %.1f", tt_hdr, ci_str, dic_i)
      else
        sprintf("%s\n%s", tt_hdr, ci_str)
      parts[length(parts) + 1L] <- paste0(
        sprintf('<g class="ma-pooled %s" data-tt="%s">', klass, esc_attr(tt)),
        sprintf('<text class="ma-rowlabel ma-pooled-label" x="%g" y="%g">%s</text>',
                LEFT_MARGIN - 10, yc + 4, esc(kind_lbl)),
        sprintf('<rect class="ma-rowhit" x="%g" y="%g" width="%g" height="%g"/>',
                PLOT_LEFT, yc - POOLED_H / 2 + 1, PLOT_W, POOLED_H - 2),
        sprintf('<polygon class="ma-diamond" points="%s"/>', pts),
        sprintf('<text class="ma-citext" x="%g" y="%g">%s</text>',
                width - RIGHT_MARGIN + 10, yc + 4, esc(ci_str)),
        '</g>'
      )
      y_run <- y_run + POOLED_H
    }
  }

  axis_y <- plot_bottom + 6
  parts[length(parts) + 1L] <- sprintf(
    '<line class="ma-axis" x1="%g" y1="%g" x2="%g" y2="%g"/>',
    PLOT_LEFT, axis_y, PLOT_LEFT + PLOT_W, axis_y)

  for (t in ticks) {
    tx <- xfn(t)
    parts[length(parts) + 1L] <- sprintf(
      '<line class="ma-tick" x1="%g" y1="%g" x2="%g" y2="%g"/>',
      tx, axis_y, tx, axis_y + 4)
    parts[length(parts) + 1L] <- sprintf(
      '<text class="ma-ticklabel" x="%g" y="%g">%s</text>',
      tx, axis_y + 16, esc(forest_tick_label(scale, t)))
  }

  if (AXIS_LABEL_H > 0) {
    if (length(axis_label) && nzchar(axis_label)) {
      parts[length(parts) + 1L] <- sprintf(
        '<text class="ma-axislabel" x="%g" y="%g">%s</text>',
        PLOT_LEFT + PLOT_W / 2, axis_y + 30, esc(axis_label))
    }
    dir_y <- axis_y + 46
    if (length(dir_left) && nzchar(dir_left)) {
      parts[length(parts) + 1L] <- sprintf(
        '<text class="ma-dir ma-dir-left" x="%g" y="%g">%s</text>',
        PLOT_LEFT, dir_y, esc(dir_left))
    }
    if (length(dir_right) && nzchar(dir_right)) {
      parts[length(parts) + 1L] <- sprintf(
        '<text class="ma-dir ma-dir-right" x="%g" y="%g">%s</text>',
        PLOT_LEFT + PLOT_W, dir_y, esc(dir_right))
    }
  }

  parts[length(parts) + 1L] <- "</svg>"
  paste(parts, collapse = "")
}
