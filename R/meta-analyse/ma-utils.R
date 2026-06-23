nma_results <- function(m, base_dist=NA, label=NA, t=NA, reft=NA) {
  
  results <- list()
  
  if (any(class(m) == "stan_nma")) {
    if(m$likelihood == "ordered") {
    
      # Generate MCMC trace for response rates
      rates <- predict(
          m, type = "response",
          baseline = base_dist,
          baseline_type = "response",
          summary = FALSE
        )$sims |> 
        posterior::as_draws_df() |>
        # Convert to long format
        pivot_longer(!starts_with("."), names_to = "param", values_to = "trace") |>
        # Extract drug name
        mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=,)"),
               endpoint = str_extract(param, pattern = "(?<=\\, ).*.(?=])")) |> 
        suppressWarnings()
      results[[1]] <- summarise(
        .by = c(drug, endpoint),
        rates,
        mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
      ) |> 
        mutate(
          type = "network",
          method = "multinomial",
          effects = m$trt_effects,
          comp_tx = drug,
          ref_tx = NA,
          measure = "rate",
        ) |> 
        select(-drug)
      
      for (i in 1:nrow(comparisons)) {
        if (!all(comparisons[i,] %in% rates$drug)) next
        pairwise <- rates |> 
          select(-param) |> 
          filter(drug %in% comparisons[i,]) |> 
          pivot_wider(names_from = drug, values_from = trace)
        pairwise$rd <- pairwise[[comparisons[i,1]]] - pairwise[[comparisons[i,2]]]
        results[[i + 1]] <- summarise(
          .by = "endpoint",
          pairwise,
          mean = mean(rd), lower = quantile(rd, 0.025), upper = quantile(rd, 0.975)
        ) |> mutate(
          type = "network",
          method = "multinomial",
          effects = m$trt_effects,
          ref_tx = comparisons[i,2],
          comp_tx = comparisons[i,1],
          measure = "rd"
        )
      }
    } else if(m$likelihood == "normal") {
      rates <- predict(
        m, type = "response",
        baseline = base_dist,
        summary = FALSE,
        baseline_type = "response"
      )$sims |> 
        posterior::as_draws_df() |>
        # Convert to long format
        pivot_longer(!starts_with("."), names_to = "param", values_to = "trace") |>
        # Extract drug name
        mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=])")) |> 
        suppressWarnings()
      results[[1]] <- summarise(
        .by = drug,
        rates,
        mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
      ) |> 
        mutate(
          type = "network",
          method = "normal",
          endpoint = label,
          effects = m$trt_effects,
          ref_tx = NA,
          comp_tx = drug,
          measure = "cfb", # Change from baseline
        ) |> 
        select(-drug)
      
      for (i in 1:nrow(comparisons)) {
        if (!all(comparisons[i,] %in% rates$drug)) next
        pairwise <- rates |> 
          select(-param) |> 
          filter(drug %in% comparisons[i,]) |> 
          pivot_wider(names_from = drug, values_from = trace)
        pairwise$rd <- pairwise[[comparisons[i,1]]] - pairwise[[comparisons[i,2]]]
        results[[i + 1]] <- summarise(
          pairwise,
          mean = mean(rd), lower = quantile(rd, 0.025), upper = quantile(rd, 0.975)
        ) |> mutate(
          type = "network",
          method = "normal",
          endpoint = label,
          effects = m$trt_effects,
          ref_tx = comparisons[i,2],
          comp_tx = comparisons[i,1],
          measure = "diff_cfb"
        )
      }
    } else if(m$likelihood == "binomial") {
      
      # Generate MCMC trace for response rates
      rates <- predict(
        m, type = "response",
        baseline = base_dist,
        summary = FALSE,
        baseline_type = "response"
      )$sims |> 
        posterior::as_draws_df() |>
        # Convert to long format
        pivot_longer(!starts_with("."), names_to = "param", values_to = "trace") |>
        # Extract drug name
        mutate(drug = str_extract(param, pattern = "(?<=\\[).*.(?=])")) |> 
        suppressWarnings()
      results[[1]] <- summarise(
        .by = drug,
        rates,
        mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
      ) |> 
        mutate(
          type = "network",
          method = "binomial",
          endpoint = label,
          effects = m$trt_effects,
          comp_tx = drug,
          ref_tx = NA,
          measure = "rate",
        ) |> 
        select(-drug)
      
      for (i in 1:nrow(comparisons)) {
        if (!all(comparisons[i,] %in% rates$drug)) next
        pairwise <- rates |> 
          select(-param) |> 
          filter(drug %in% comparisons[i,]) |> 
          pivot_wider(names_from = drug, values_from = trace)
        pairwise$rd <- pairwise[[comparisons[i,1]]] - pairwise[[comparisons[i,2]]]
        results[[i + 1]] <- summarise(
          pairwise,
          mean = mean(rd), lower = quantile(rd, 0.025), upper = quantile(rd, 0.975)
        ) |> mutate(
          type = "network",
          method = "binomial",
          endpoint = label,
          effects = m$trt_effects,
          ref_tx = comparisons[i,2],
          comp_tx = comparisons[i,1],
          measure = "rd"
        )
      }
    }
  } else if(any(class(m) == "metaprop")) {
    results[[1]] <- data.frame(
      endpoint = label,
      type = "univariate",
      method = "logit",
      effects = "fixed",
      ref_tx = NA,
      comp_tx = t,
      measure = "rate",
      mean = plogis(m$TE.fixed),
      lower = plogis(m$lower.fixed),
      upper = plogis(m$upper.fixed)
    )
    results[[2]] <- data.frame(
      endpoint = label,
      type = "univariate",
      method = "logit",
      effects = "random",
      ref_tx = NA,
      comp_tx = t,
      measure = "rate",
      mean = plogis(m$TE.random),
      lower = plogis(m$lower.random),
      upper = plogis(m$upper.random)
    )
  } else if(any(class(m) == "metabin")) {
    results[[1]] <- data.frame(
      endpoint = label,
      type = "pairwise",
      method = "binomial",
      effects = "fixed",
      ref_tx = reft,
      comp_tx = t,
      measure = "rd",
      mean = m$TE.fixed,
      lower = m$lower.fixed,
      upper = m$upper.fixed
    )
    results[[2]] <- data.frame(
      endpoint = label,
      type = "pairwise",
      method = "binomial",
      effects = "random",
      ref_tx = reft,
      comp_tx = t,
      measure = "rd",
      mean = m$TE.random,
      lower = m$lower.random,
      upper = m$upper.random
    )
  } else if (any(class(m) == "metacont")) {
    results[[1]] <- data.frame(
      endpoint = label,
      type = "pairwise",
      method = "normal",
      effects = "fixed",
      ref_tx = reft,
      comp_tx = t,
      measure = "diff_cfb",
      mean = m$TE.fixed,
      lower = m$lower.fixed,
      upper = m$upper.fixed
    )
    results[[2]] <- data.frame(
      endpoint = label,
      type = "pairwise",
      method = "normal",
      effects = "random",
      ref_tx = reft,
      comp_tx = t,
      measure = "diff_cfb",
      mean = m$TE.random,
      lower = m$lower.random,
      upper = m$upper.random
    )
  } else if (any(class(m) == "metagen")) {
    results[[1]] <- data.frame(
      endpoint = label,
      type = "univariate",
      method = "normal",
      effects = "fixed",
      ref_tx = NA,
      comp_tx = t,
      measure = "cfb",
      mean = m$TE.fixed,
      lower = m$lower.fixed,
      upper = m$upper.fixed
    )
    results[[2]] <- data.frame(
      endpoint = label,
      type = "univariate",
      method = "normal",
      effects = "random",
      ref_tx = NA,
      comp_tx = t,
      measure = "cfb",
      mean = m$TE.random,
      lower = m$lower.random,
      upper = m$upper.random
    )
  }
  
  bind_rows(results)
}

beta_dist_metaprop <- function(mod, effects) {
  
  if (effects == "fixed") {
    mu = plogis(mod$TE.fixed)
    lower = plogis(mod$lower.fixed)
    upper = plogis(mod$upper.fixed)
    se = (upper - lower) / (2 * 1.96)
  } else if (effects == "random") {
    mu = plogis(mod$TE.random)
    lower = plogis(mod$lower.random)
    upper = plogis(mod$upper.random)
    se = (upper - lower) / (2 * 1.96)
  }
  
  var <- se^2
  
  # Constraint check
  if (var >= mu * (1 - mu)) {
    stop("Variance (SE^2) is too high for a valid Beta distribution.")
  }
  
  # Calculate the common factor
  factor <- (mu * (1 - mu) / var) - 1
  
  alpha <- mu * factor
  beta <- (1 - mu) * factor
  
  return(distr(qbeta, alpha, beta))
}

