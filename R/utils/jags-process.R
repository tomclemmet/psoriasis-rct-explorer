library(dplyr)
library(stringr)


process_jags <- function(mod) {
  if (inherits(mod, "jags_nma_fit")) return(mod)
  
  param_lookup <- data.frame(
    label = pasi_id_lookup$drug,
    param = paste0("d[", seq(1:length(pasi_id_lookup$drug)), "]")
  ) |> bind_rows(distinct(pasi_id_lookup, class, cl) |> 
                mutate(label = class, param = paste0("m[", cl, "]")) |>
                select(-class, -cl))
  
  fitted_values <- mod$BUGSoutput$summary |>
    as_tibble(rownames = "param") |> 
    select(param, mean, `2.5%`, `97.5%`, Rhat) |> 
    filter(str_starts(param, "rhat\\[")) |> 
    separate_wider_regex(
      param, 
      patterns = c(".*\\[", id_row = "\\d+", ",\\s*", id_arm = "\\d+", ",\\s*", 
                   id_cat = "\\d+", "\\].*")
    ) |> 
    mutate(across(starts_with("id_"), as.integer)) |> 
    rename(fitted = mean) |> 
    mutate(actual = pasi_jags$r[cbind(id_row, id_arm, id_cat)]) |> 
    select(starts_with("id"), fitted, actual)
  
  out <- list(
    results = mod$BUGSoutput$summary |>
      as_tibble(rownames = "param") |> 
      left_join(param_lookup, by = "param") |> 
      relocate(label, .after = param) |> 
      as.data.frame(),
    
    trace = posterior::as_draws_df(mod$BUGSoutput$sims.array) |> 
      select(starts_with("."), starts_with(c("d[", "z[", "sd", "B", "mu", "totresdev"))),
    
    dev_table = mod$BUGSoutput$summary |>
      as_tibble(rownames = "param") |> 
      as.data.frame() |> 
      select(param, mean, `2.5%`, `97.5%`, Rhat) |> 
      filter(str_starts(param, "dv\\[")) |> 
      separate_wider_regex(
        param, 
        patterns = c(".*\\[", id_row = "\\d+", ",\\s*", id_arm = "\\d+", ",\\s*", 
                     id_cat = "\\d+", "\\].*")
      ) |> 
      mutate(across(starts_with("id_"), as.integer)) |> 
      left_join(row_lookup, by = c("id_row" = "row", "id_arm" = "arm_no"),
                relationship = "many-to-one") |> 
      arrange(ref_id, id_arm, id_cat) |> 
      left_join(fitted_values, by = c("id_row", "id_arm", "id_cat")) |> 
      relocate(ref_id, trial, drug, id_cat, fitted, actual, mean, `2.5%`, `97.5%`,
               Rhat, id_row, id_arm, t),
    
    totresdev = mod$BUGSoutput$mean$totresdev,
    
    pV = mod$BUGSoutput$pV,
    
    DIC = as.numeric(mod$BUGSoutput$mean$totresdev + mod$BUGSoutput$pV)
  )
  
  out$summary <- out$results |> 
    filter(str_starts(param, "d\\[|sd|m\\[|z|B_|mubar"), !str_detect(param, ",")) |> 
    arrange(str_detect(param, "B|mubar"))
  
  class(out) <- c("jags_nma_fit", class(out))
  
  out
}

print.jags_nma_fit <- function(m) {
  totresdev <- m$results[m$results$param == "totresdev", 3]
  
  message(paste0("totresdev = ", round(totresdev, 1), " on ", nrow(m$dev_table), " data points, pV = ", round(m$pV, 1), ", DIC = ", round(m$DIC, 1)))
  print(m$summary)
  invisible(m$summary)
  
}

# Function to display goodness-of-fit stats for a list of models
check_fits <- function(mods) {
  for (i in 1:length(mods)) {
    m <- process_jags(mods[[i]])
    totresdev <- m$results[m$results$param == "totresdev", 3]
    message(paste0(names(mods)[i], ": totresdev = ", round(totresdev, 1), " on ", nrow(m$dev_table), " data points, pV = ", round(m$pV, 1), ", DIC = ", round(m$DIC, 1)))
  }
}

# Function to compare model outputs given a list of jags models
compare_jags <- function(mods) {
  
  n_params <- c()
  
  for (i in 1:length(mods)) {
    n_params[i] <- nrow(process_jags(mods[[i]])$summary)
  }
  
  tab <- data.frame(
    param = process_jags(mods[[which.max(n_params)]])$summary$param,
    drug = process_jags(mods[[which.max(n_params)]])$summary$label
  )
  
  for (i in 1:length(mods)) {
    coefs <- process_jags(mods[[i]])$summary |> 
      mutate(meansd = paste0(round(mean, 3), " (", round(sd, 3), "; ", round(Rhat, 3), ")")) |> 
      select(param, meansd)
      
    
    tab <- tab |> left_join(coefs, by = "param")
  }
  
  totresdev <- lapply(mods, \(x) {as.character(round(process_jags(x)$totresdev, 3))}) |> 
    as.data.frame() |> 
    mutate(param = "totresdev")
  pV <- lapply(mods, \(x) {as.character(round(process_jags(x)$pV, 3))}) |> 
    as.data.frame() |> 
    mutate(param = "pV")
  dic <- lapply(mods, \(x) {as.character(round(process_jags(x)$DIC, 3))}) |> 
    as.data.frame() |> 
    mutate(param = "dic")

  
  tab |> 
    rename_with(~ names(mods), .cols = starts_with("meansd")) |>
    bind_rows(totresdev, pV, dic) |> 
    arrange(desc(param %in% c("totresdev", "pV", "dic")))
}
