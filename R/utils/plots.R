devplot <- function(m1, m2, xlab = "Model 1", ylab = "Model 2", output = c("plot", "table")) {
  output = match.arg(output)
  
  devdev <- inner_join(process_jags(m1)$dev_table, process_jags(m2)$dev_table, 
                       by = c("id_row", "id_arm", "id_cat", "trial",
                              "ref_id", "t", "drug")) |> 
    mutate(.by = ref_id, diff = mean.x - mean.y, .after = id_cat)
  
  if (output == "plot") {
    ggplot(devdev, aes(x = mean.x, y = mean.y)) +
      geom_point(alpha = 0.5, shape = 16) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "blue") +
      geom_abline(intercept = -0.5, slope = 1, linetype = 3, colour = "red") +
      geom_abline(intercept = -1, slope = 1, linetype = 3, colour = "red") +
      theme_classic() +
      labs(title = "Deviance-deviance plot", x = xlab, y = ylab) +
      scale_color_viridis_d()
  } else if (output == "table") {
    devdev
  }
  
}

forest <- function(m) {
  df <- m$summary |> 
    filter(label %notin% c("Mirikizumab", "Phototherapy", "Xeligekimab", "Netakimab", "Roflumilast", "Icotrokinra", "Tofacitinib"), param != "mubar") |> 
    mutate(group = substr(param, 1, 1), group = factor(if_else(group == "s", "sd", group), levels = c("d", "z", "sd", "B")), label = if_else(is.na(label), param, label), rank = if_else(group == "d", mean, NA)) |> 
    group_by(group) |> 
    arrange(desc(rank), .by_group = TRUE)
  df$label <- factor(df$label, levels = rev(df$label))
  ggplot(df, aes(y = label, x = mean)) +
    geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`), shape = 15) +
    facet_grid(group ~ ., scales = "free", space = "free") +
    theme_bw()
}

forests <- function(..., lab = NA, prob = FALSE) {
  mods <- list(...)
  nms <- names(mods)
  
  if (is.null(nms)) {
    exprs <- as.list(substitute(list(...)))[-1]
    nms <- vapply(exprs, deparse, character(1))
    nms <- sub("^.*\\$", "", nms)   # strip everything up to and including "$"
  }
  names(mods) <- nms
  
  if (length(mods) == 1 && is.list(mods[[1]])) {
    mods <- mods[[1]]
  }
  df <- list()
  
  if (!prob) {
    
    for (i in 1:length(mods)) {
      df[[names(mods)[i]]] <- mods[[i]]$summary |>
        filter(label %notin% c(
          "Icotrokinra", "Mirikizumab", "Netakimab", "Orismilast", "Roflumilast",
          "Phototherapy", "Sonelokimab", "Tofacitinib", "Xeligekimab", "Zasocitinib"
        ), param != "mubar") |>
        mutate(group = substr(param, 1, 1), group = factor(if_else(group == "s", "sd", group), levels = c("d", "z", "sd", "B")), label = if_else(is.na(label), param, label), rank = if_else(group == "d", mean, NA)) |>
        group_by(group) |>
        arrange(desc(rank), .by_group = TRUE)
      df[[names(mods)[i]]]$label <- factor(df[[names(mods)[i]]]$label, levels = rev(df[[names(mods)[i]]]$label))
    }
    
    if (!is.na(lab)) {
      df <- bind_rows(df, .id = "model") |> filter(label %in% lab)
    } else {
      df <- bind_rows(df, .id = "model")
    }
    
    ggplot(df, aes(y = label, x = mean)) +
      geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = model), shape = 15, size = 0.2, position = position_dodge(0.5)) +
      facet_grid(group ~ ., scales = "free", space = "free") +
      theme_bw()
    
  } else {
    outcomes <- c("pasi50", "pasi75", "pasi90", "pasi100")
    
    for (i in 1:length(mods)) {
      drug_rank <- mods[[i]]$results |> 
        filter(str_starts(param, "prob")) |> 
        mutate(drug = pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))]) |> 
        slice_head(n = 1, by = drug) |> 
        arrange(mean)
      df[[names(mods)[i]]] <- mods[[i]]$results |>
        filter(str_starts(param, "prob")) |> 
        mutate(
          drug = factor(
            pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))],
            levels = drug_rank$drug
          ),
          outcome = factor(
            outcomes[as.numeric(str_extract(param, "(?<=\\[).*?(?=,)"))],
            levels = c("pasi50", "pasi75", "pasi90", "pasi100")
          )
        ) |> 
        filter(drug %notin% c(
          "Icotrokinra", "Mirikizumab", "Netakimab", "Orismilast", "Roflumilast",
          "Phototherapy", "Sonelokimab", "Tofacitinib", "Xeligekimab", "Zasocitinib"
        ))
      df[[names(mods)[i]]]$label <- factor(df[[names(mods)[i]]]$label, levels = rev(df[[names(mods)[i]]]$label))
    }
    
    if (!is.na(lab)) {
      df <- bind_rows(df, .id = "model") |> filter(drug %in% lab)
    } else {
      df <- bind_rows(df, .id = "model")
    }
    
    ggplot(df, aes(y = drug, x = mean)) +
      geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`, colour = outcome, shape = model), size = 0.2, position = position_dodge(0.5)) +
      theme_bw()
    
  }
}

converge <- function(m, param, xmin = 0) {
  
  m$trace |> 
    select(starts_with("."), all_of(param)) |> 
    group_by(.chain) |> 
    mutate(cummean = cumsum(.data[[param]]) / .iteration) |> 
    filter(.iteration >= xmin) |> 
    ggplot() +
    geom_line(aes(x = .iteration, y = cummean, colour = factor(.chain))) +
    theme_bw()
  
}

