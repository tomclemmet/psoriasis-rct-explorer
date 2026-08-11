library(stringr)

nma_results <- function(m, base_dist=NA, method = "standard", effects = NA, label=NA, t=NA, reft=NA) {
  
  results <- list()
  drug_order <- c(
      "Placebo", "Acitretin", "Adalimumab", "Apremilast", "Bimekizumab",
      "Brodalumab", "Certolizumab", "Cyclosporin", "Deucravacitinib", "Etanercept",
      "Fumaric acid esters", "Guselkumab", "Icotrokinra", "Infliximab", "Ixekizumab",
      "Izokibep", "Methotrexate", "Mirikizumab", "Netakimab", "Orismilast",
      # "Phototherapy", 
      "Risankizumab", "Roflumilast", "Secukinumab", "Sonelokimab",
      "Tildrakizumab", "Tofacitinib", "Ustekinumab", "Xeligekimab", "Zasocitinib"
    )
  thresholds <- c("pasi50", "pasi75", "pasi90", "pasi100")
  
  if (any(class(m) == "rjags")) {
    ##!!!!! CURRENTLY IGNORES BASE DIST
    dic <- process_jags(m)$DIC
    
    # Generate MCMC traces for response rates
    rates <- process_jags(m)$trace |> 
      select(starts_with(c("prob", "."))) |> 
      # Convert to long format
      pivot_longer(starts_with("prob"), names_to = "param", values_to = "trace") |> 
      # Extract drug name
      mutate(drug = drug_order[as.numeric(str_extract(param, pattern = "(?<=,\\s?)\\d+(?=\\])"))],
             endpoint = thresholds[as.numeric(str_extract(param, pattern = "(?<=prob\\[)\\d+(?=,)"))]) |> 
      suppressWarnings() 
    
    # Summarise response rates for each treatment
    results[[1]] <- summarise(
      .by = c(drug, endpoint),
      rates,
      mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
    ) |> 
      mutate(
        type = "network",
        likelihood = "multinomial",
        method = method,
        effects = effects,
        comp_tx = drug,
        ref_tx = NA,
        measure = "rate",
        dic = dic
      ) |> 
      select(-drug)
    
    # Summarise risk differences for each comparison
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
        likelihood = "multinomial",
        method = method,
        effects = effects,
        ref_tx = comparisons[i,2],
        comp_tx = comparisons[i,1],
        measure = "rd",
        dic = dic
      )
    }
    
  } else {
  
    if (any(class(m) == "stan_nma")) {
      dic <- dic(m)$dic
      
      if(m$likelihood == "ordered") {
      
        # Generate MCMC traces for response rates
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
        
        # Summarise response rates for each treatment
        results[[1]] <- summarise(
          .by = c(drug, endpoint),
          rates,
          mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
        ) |> 
          mutate(
            type = "network",
            likelihood = "multinomial",
            method = method,
            effects = m$trt_effects,
            comp_tx = drug,
            ref_tx = NA,
            measure = "rate",
            dic = dic
          ) |> 
          select(-drug)
        
        # Summarise risk differences for each comparison
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
            likelihood = "multinomial",
            method = method,
            effects = m$trt_effects,
            ref_tx = comparisons[i,2],
            comp_tx = comparisons[i,1],
            measure = "rd",
            dic = dic
          )
        }
      
      } else if(m$likelihood == "normal") {
        
        # Generate MCMC traces for response rates
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
        
        # Summarise response rates for each treatment
        results[[1]] <- summarise(
          .by = drug,
          rates,
          mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
        ) |> 
          mutate(
            type = "network",
            likelihood = "normal",
            method = method,
            endpoint = label,
            effects = m$trt_effects,
            ref_tx = NA,
            comp_tx = drug,
            measure = "cfb", # Change from baseline
            dic = dic
          ) |> 
          select(-drug)
        
        # Summarise risk differences for each comparison
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
            likelihood = "normal",
            method = method,
            endpoint = label,
            effects = m$trt_effects,
            ref_tx = comparisons[i,2],
            comp_tx = comparisons[i,1],
            measure = "diff_cfb",
            dic = dic
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
        
        # Summarise response rates for each treatment
        results[[1]] <- summarise(
          .by = drug,
          rates,
          mean = mean(trace), lower = quantile(trace, 0.025), upper = quantile(trace, 0.975)
        ) |> 
          mutate(
            type = "network",
            likelihood = "binomial",
            method = method,
            endpoint = label,
            effects = m$trt_effects,
            comp_tx = drug,
            ref_tx = NA,
            measure = "rate",
            dic = dic
          ) |> 
          select(-drug)
        
        # Summarise risk differences for each comparison
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
            likelihood = "binomial",
            method = method,
            endpoint = label,
            effects = m$trt_effects,
            ref_tx = comparisons[i,2],
            comp_tx = comparisons[i,1],
            measure = "rd",
            dic = dic
          )
        }
      }
      
    # Process pairwise or univariate meta-analysis results
    } else if(any(class(m) == "metaprop")) {
      results[[1]] <- data.frame(
        endpoint = label,
        type = "univariate",
        likelihood = "logit",
        method = method,
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
        likelihood = "logit",
        method = method,
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
        likelihood = "binomial",
        method = method,
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
        likelihood = "binomial",
        method = method,
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
        likelihood = "normal",
        method = method,
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
        likelihood = "normal",
        method = method,
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
        likelihood = "normal",
        method = method,
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
        likelihood = "normal",
        method = method,
        effects = "random",
        ref_tx = NA,
        comp_tx = t,
        measure = "cfb",
        mean = m$TE.random,
        lower = m$lower.random,
        upper = m$upper.random
      )
    }
  }
  bind_rows(results)
}

process_jags <- function(mod) {
  drug_lookup <- data.frame(
    drug = c(
      "Placebo", "Acitretin", "Adalimumab", "Apremilast", "Bimekizumab",
      "Brodalumab", "Certolizumab", "Cyclosporin", "Deucravacitinib", "Etanercept",
      "Fumaric acid esters", "Guselkumab", "Icotrokinra", "Infliximab", "Ixekizumab",
      "Izokibep", "Methotrexate", "Mirikizumab", "Netakimab", "Orismilast",
      # "Phototherapy",
      "Risankizumab", "Roflumilast", "Secukinumab", "Sonelokimab",
      "Tildrakizumab", "Tofacitinib", "Ustekinumab", "Xeligekimab", "Zasocitinib"
    ),
    index = paste0("d[", seq(1:length(drug_order)), "]")
  )
  
  list(
    summary = mod$BUGSoutput$summary |> 
      as_tibble(rownames = "param") |> 
      left_join(drug_lookup, by = c("param" = "index")) |> 
      relocate(drug, .after = param) |> 
      filter(!str_detect(param, "prob")) |> 
      arrange(param %in% c("deviance", "totresdev")) |> 
      as.data.frame(),
    trace = posterior::as_draws_df(mod$BUGSoutput$sims.array),
    totresdev = mod$BUGSoutput$mean$totresdev,
    pV = mod$BUGSoutput$pV,
    DIC = as.numeric(mod$BUGSoutput$mean$totresdev + mod$BUGSoutput$pV)
  )
}

# Function to compare model outputs given a list of jags models
compare_jags <- function(mods) {
  
  n_params <- c()
  
  for (i in 1:length(mods)) {
    n_params[i] <- nrow(process_jags(mods[[i]])$summary)
  }
  
  tab <- data.frame(
    param = process_jags(mods[[which.max(n_params)]])$summary$param,
    drug = process_jags(mods[[which.max(n_params)]])$summary$drug
  )
  
  for (i in 1:length(mods)) {
    coefs <- process_jags(mods[[i]])$summary |> 
      mutate(meansd = paste0(round(mean, 3), " (", round(sd, 3), ")")) |> 
      select(param, meansd)
      
    
    tab <- tab |> left_join(coefs, by = "param")
  }
  
  dic <- lapply(mods, \(x) {as.character(round(process_jags(x)$DIC, 3))}) |> 
    as.data.frame() |> 
    mutate(param = "dic")
  pV <- lapply(mods, \(x) {as.character(round(process_jags(x)$pV, 3))}) |> 
    as.data.frame() |> 
    mutate(param = "pV")
  
  tab |> 
    rename_with(~ names(mods), .cols = starts_with("meansd")) |>
    bind_rows(dic, pV)
}

# Helper function to produce a beta multinma::distr() object to use as a 
# baseline in multinma::predict() based on a fitted meta::metaprop() object
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



