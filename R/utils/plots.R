# Functions to produce plots for processed NMA results
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
source("R/utils/jags-process.R")


# Deviance comparison plot for 2 models, can also output a table
devplot <- function(m1, m2, xlab = "Model 1", ylab = "Model 2", output = c("plot", "table")) {
  output = match.arg(output)
  
  # Merge deviance tables from each model and calculate difference
  devdev <- inner_join(process_jags(m1)$dev_table, process_jags(m2)$dev_table, 
                       by = c("id_row", "id_arm", "id_cat", "trial",
                              "ref_id", "t", "drug")) |> 
    mutate(.by = ref_id, diff = mean.x - mean.y, .after = id_cat)
  
  if (output == "plot") {
    # Produce plot with 45-degree guide lines
    ggplot(devdev, aes(x = mean.x, y = mean.y)) +
      geom_point(alpha = 0.5, shape = 16) +
      geom_text_repel(data = slice_max(devdev, diff, n = 5), 
                      aes(label = paste0(trial, " [", ref_id, "]")),
                      nudge_x = 0.15, nudge_y = -0.15,
                      hjust = 0, vjust = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "blue") +
      geom_abline(intercept = -0.5, slope = 1, linetype = 3, colour = "red") +
      geom_abline(intercept = -1, slope = 1, linetype = 3, colour = "red") +
      theme_classic() +
      labs(title = "Deviance-deviance plot", x = xlab, y = ylab) +
      scale_color_viridis_d()
  } else if (output == "table") {
    # If desired, simply output the table
    devdev
  }
  
}

# Simple forest plot for key parameters
forest <- function(m) {
  
  m <- process_jags(m)
  
  # Filter unwanted drugs, then group and order parameters for plotting
  df <- m$summary |> 
    filter(label %notin% c("Mirikizumab", "Phototherapy", "Xeligekimab", "Netakimab", "Roflumilast", "Icotrokinra", "Tofacitinib", "Izokibep"), param != "mubar") |> 
    mutate(group = substr(param, 1, 1), group = factor(if_else(group == "s", "sd", group), levels = c("d", "z", "sd", "B")), label = if_else(is.na(label), param, label), rank = if_else(group == "d", mean, NA)) |> 
    group_by(group) |> 
    arrange(desc(rank), .by_group = TRUE)
  df$label <- factor(df$label, levels = rev(df$label))
  
  # Produce plot
  ggplot(df, aes(y = label, x = mean)) +
    geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`), shape = 15) +
    facet_grid(group ~ ., scales = "free", space = "free") +
    theme_bw()
}

# Forest plot for multiple models, option to display predicted probabilities or 
# parameter values and filter for specific parameters
forests <- function(..., lab = NA, grp = NA, remove = NA) {
  mods <- list(...)
  
  # Handle case where a single list is passed
  if (length(mods) == 1 && !inherits(mods[[1]], "jags_nma_fit")) {
    mods <- mods[[1]]
  }
  
  nms <- names(mods)
  
  # Assign variable names as list labels if passed as separate items
  if (is.null(nms)) {
    exprs <- as.list(substitute(list(...)))[-1]
    nms <- vapply(exprs, deparse, character(1))
    nms <- sub("^.*\\$", "", nms)
  }
  names(mods) <- nms
  
  mods <- lapply(mods, process_jags)
  
  df <- list()
  
  for (i in 1:length(mods)) {
    
    # For each model, filter unwanted drugs, then group and order parameters
    df[[names(mods)[i]]] <- mods[[i]]$summary |>
      filter(label %notin% c(
        "Icotrokinra", "Mirikizumab", "Netakimab", "Orismilast", "Roflumilast",
        "Phototherapy", "Sonelokimab", "Tofacitinib", "Xeligekimab", "Zasocitinib"
      ), param != "mubar") |>
      mutate(
        group = factor(case_when(
          str_detect(param, "d\\[") ~ "Treatment effect",
          str_detect(param, "m\\[") ~ "Class effect",
          str_detect(param, "sd") ~ "Heterogeneity",
          str_detect(param, "z\\[") ~ "Multinomial cutpoint",
          str_detect(param, "B|beta") ~ "Regression coefficient"
          ), level = c("Treatment effect", "Class effect", "Heterogeneity", 
                       "Multinomial cutpoint", "Regression coefficient")),
        label = if_else(is.na(label), param, label),
        rank = if_else(group %in% c("Treatment effect", "Class effect"), mean, NA),
        converged = Rhat < 1.1,
        `2.5%` = if_else(!converged, NA, `2.5%`),
        `97.5%` = if_else(!converged, NA, `97.5%`)
      ) |>
      group_by(group) |>
      arrange(desc(rank), .by_group = TRUE)
    df[[names(mods)[i]]]$label <- factor(df[[names(mods)[i]]]$label, levels = rev(df[[names(mods)[i]]]$label))
  }
  
  # Select chosen parameter if specified
  if (any(!is.na(lab))) {
    df <- bind_rows(df, .id = "model") |> filter(label %in% lab)
  } else if (any(!is.na(remove))) {
    df <- bind_rows(df, .id = "model") |> filter(!(label %in% remove))
  }
    df <- bind_rows(df, .id = "model")
  }
  
  # Define model ordering for plot
  if (any(df$model == "unadjusted")) {
    df <- df |> mutate(model = factor(model, levels = c("unadjusted", sort(setdiff(unique(model), "unadjusted")))))
  }
    
  # Select chosen group(s) if specified
  if (any(!is.na(grp))) {
    df <- df |> filter(group %in% grp)
  }
  
  # Produce plot
  ggplot(df, aes(y = label, x = mean)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = model, shape = converged), size = 0.2, position = position_dodge(0.5, reverse = TRUE)) +
    facet_wrap(~ group, scales = "free", space = "free_y") +
    theme_bw() +
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 4))
  
}

# Stacked bar chart for different leveks of PASI response
props_plot <- function(m) {
  
  outcomes <- c(
    "pasi50", "pasi75", "pasi90", "pasi100"
  )
  
  # Rank by predicted PASI 50 rate
  drug_rank <- m$results |> 
    filter(str_starts(param, "prob")) |> 
    mutate(drug = pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))]) |> 
    slice_head(n = 1, by = drug) |> 
    arrange(mean)
  
  # Extract predicted probabilities, then label and rank accordingly and produce
  # plot
  m$results |> 
    filter(str_starts(param, "prob")) |> 
    mutate(
      drug = factor(
        pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))],
        levels = drug_rank$drug
      ),
      outcome = factor(
        outcomes[as.numeric(str_extract(param, "(?<=\\[).*?(?=,)"))],
        levels = c("pasi50", "pasi75", "pasi90", "pasi100"),
        labels = c("PASI 50-75", "PASI 75-90", "PASI 90-100", "PASI 100")
      )
    ) |> 
    arrange(drug, desc(outcome)) |> 
    mutate(.by = drug, mean = mean - lag(mean, default = 0), .after = mean) |> 
    mutate(outcome = forcats::fct_rev(outcome)) |> 
    filter(drug %notin% c(
      "Icotrokinra", "Mirikizumab", "Netakimab", "Orismilast", "Roflumilast", 
      "Phototherapy", "Sonelokimab", "Tofacitinib", "Xeligekimab", "Zasocitinib")
    ) |> 
    ggplot(aes(x = mean, y = drug)) +
    geom_col(aes(fill = outcome), position = position_stack(reverse = TRUE)) +
    theme_classic() +
    scale_fill_viridis_d() +
    theme(legend.position = "right", axis.title.y = element_blank(), legend.title = element_blank(),
          plot.title = element_text(face = "bold")) +
    labs(title = "Predicted PASI Response", x = "Proportion")
}



