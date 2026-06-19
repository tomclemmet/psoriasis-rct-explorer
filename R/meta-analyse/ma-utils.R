nma_results <- function(m, base_mean=NA, base_se=NA, label=NA, t=NA, reft=NA) {
  
  results <- list()
  
  if (any(class(m) == "stan_nma")) {
    if(m$likelihood == "ordered") {
    
      # Generate MCMC trace for response rates
      rates <- predict(
          m, type = "response",
          baseline = distr(qnorm, base_mean, base_se),
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
          effects = m$trt_effects,
          ref_tx = comparisons[i,2],
          comp_tx = comparisons[i,1],
          measure = "rd"
        )
      }
    } else if(m$likelihood == "normal") {
      rates <- predict(
        m, type = "response",
        baseline = distr(qnorm, base_mean, base_se),
        summary = FALSE
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
          endpoint = label,
          effects = m$trt_effects,
          ref_tx = comparisons[i,2],
          comp_tx = comparisons[i,1],
          measure = "diff_cfb"
        )
      }
    }
  } else if(any(class(m) == "metaprop")) {
    results[[1]] <- data.frame(
      endpoint = label,
      type = "univariate",
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
