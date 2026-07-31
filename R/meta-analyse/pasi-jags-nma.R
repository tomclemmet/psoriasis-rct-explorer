library(R2jags)

pso_jags <- function(
    data, 
    filename = "JAGS/temp.jags",
    effects = c("fixed", "random"),
    cutpoints = c("fixed", "random"), 
    baseline = c("unadjusted", "adjusted")
) {
  effects <- match.arg(effects)
  cutpoints <- match.arg(cutpoints)
  baseline <- match.arg(baseline)
  
  if (baseline == "unadjusted") {
    data$mmu <- NULL
  }
  
  setup <- "
model {
  # *** PROGRAM STARTS
  for(i in 1:ns){                                                               # LOOP THROUGH STUDIES
    w[i, 1] <- 0                                                                # adjustment for multi-arm trials is zero for control arm
    delta[i, 1] <- 0                                                            # treatment effect is zero for control arm
    mu[i] ~ dnorm(0, .0001)                                                     # vague priors for all trial baselines
    for (k in 1:na[i]) {                                                        # LOOP THROUGH ARMS
      p[i, k, 1] <- 1                                                           # Pr(PASI >0)
      for (j in 1:(nc[i] - 1)) {                                                # LOOP THROUGH CATEGORIES
        r[i, k, j] ~ dbin(q[i, k, j], n[i, k, j])                               # binomial likelihood
        q.raw[i, k, j] <- 1 - (p[i, k, C[i, j + 1]] /                           # conditional probabilities
        max(p[i, k, C[i, j]], 1e-10)) 
        q[i, k, j] <- max(1e-7, min(1 - 1e-7, q.raw[i, k, j]))
        "
  theta <- "theta[i, k, j] <- mu[i] - delta[i, k] + "
  
  z_1d <- "z[C[i, j + 1] - 1]"
  z_2d <- "zeta[t[i, k], C[i, j + 1] - 1]"
  
  baseline_adj <- " + beta * (mu[i] - mmu) * (1 - equals(k, 1))"
  
  deviance <- "
        rhat[i, k, j] <- q[i, k, j] * n[i, k, j]                                  # predicted number events
        dv[i, k, j] <- 2 * (                                                    # Deviance contribution of each category  
          r[i, k, j] * 
            (log(max(r[i, k, j], 1e-10)) - log(max(rhat[i, k, j], 1e-10))) +            
            (n[i, k, j] - r[i, k, j]) * 
            (log(max(n[i, k, j] - r[i, k, j], 1e-10)) - 
               log(max(n[i, k, j] - rhat[i, k, j], 1e-10)))
        )
      }
      dev[i, k] <- sum(dv[i, k, 1:(nc[i] - 1)])                                 # deviance contribution of each arm"
  
  phi <- "
      for (j in 2:nc[i]) {                                                      # LOOP THROUGH CATEGORIES
        p[i, k, C[i, j]] <- 1 - phi.adj[i, k, j]                                # link function
        # adjust link function phi(x) for extreme values that can give numerical errors
        # when x< -5, phi(x)=0, when x> 5, phi(x)=1
        phi.adj[i, k, j] <- step(5 + theta[i, k, j - 1]) * 
          (step(theta[i, k, j - 1] - 5) +
             step(5 - theta[i, k, j - 1]) * phi(theta[i, k, j - 1]) )
      }
    }"
  
  delta_re <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] ~ dnorm(md[i, k], taud[i, k])
      md[i, k] <- d[t[i, k]] - d[t[i, 1]] + sw[i, k]                            # mean of LHR distributions, with multi-arm trial correction
      taud[i, k] <- tau * 2 * (k - 1) / k                                       # precision of LHR distributions (with multi-arm trial correction)
      w[i, k] <- (delta[i, k] - d[t[i, k]] + d[t[i, 1]])                        # adjustment, multi-arm RCTs
      sw[i, k] <- sum(w[i, 1:(k-1)]) / (k-1)                                    # cumulative adjustment for multi-arm trials
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                            # summed residual deviance contribution for this trial
  }"
  
  delta_fe<- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] <- d[t[i, k]] - d[t[i, 1]]
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                            # summed residual deviance contribution for this trial
  }"
  
  fez <- "
  z[1] <- 0                                                                     # set z50=0
  for (j in 2:(Cmax-1)) {                                                         # Set priors for z, for any number of categories
    z.aux[j] ~ dunif(0,5)                                                       # priors
    z[j] <- z[j-1] + z.aux[j]                                                   # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
  }"
  
  rez <- "                                                                    
  for (i in 1:ns) {zeta[i, 1] <- 0 } # set z50=0
  z[1] <- 0
  for (j in 2:(Cmax-1)) {                                                         # Set priors for z, for any number of categories
    z.aux[j] ~ dunif(0,5)                                                       # priors
    z[j] <- z[j - 1] + z.aux[j]
    for (i in 1:nt) {
      zeta.aux[i, j] ~ dnorm(z.aux[j], tauz)
      zeta[i, j] <- zeta[i, j - 1] + zeta.aux[i, j]                                              # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
    }
  }"
  
  priors <- "
  totresdev <- sum(resdev[])                                                    # Total Residual Deviance
  d[1] <- 0                                                                     # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) }                                      # vague priors for treatment effects
  beta ~ dnorm(0,.0001)
  sd ~ dunif(0, 5)                                                               # vague prior for between-trial SD
  sdz ~ dunif(0, 5)
  tau <- pow(sd,-2)   
  tauz <- pow(sdz, -2) # between-trial precision = (1 / between-trial variance)
  mubar <- mean(mu[])
  
  
  A ~ dnorm(1.097,123) 
  # calculate prob of achieving PASI 50/75/90/100 on treatment k"

  probs_fez <- "
  for (k in 1:nt) {
    for (j in 1:(Cmax - 1)) { prob[j,k] <- 1 - phi(A - d[k] + z[j]) }
  } 
  # *** PROGRAM ENDS 
}"
  
probs_rez <- "
  for (k in 1:nt) {
    for (j in 1:(Cmax - 1)) { prob[j,k] <- 1 - phi(A - d[k] + zeta[k, j]) }
  } 
  # *** PROGRAM ENDS 
}"
  
  model_code <- paste0(
    setup,
    theta,
    if (cutpoints == "fixed") z_1d else z_2d,
    if (baseline == "unadjusted") "" else baseline_adj,
    deviance,
    phi,
    if(effects == "fixed") delta_fe else delta_re,
    if (cutpoints == "fixed") fez else rez,
    priors,
    if (cutpoints == "fixed") probs_fez else probs_rez
  )
  
  writeLines(model_code, filename)
  
  params <- c(
    "d", "z", "prob",
    if (effects == "random") "sd" else NULL,
    if (cutpoints == "random") "sdz" else NULL,
    if (baseline == "adjusted") c("beta", "mubar") else NULL,
    "totresdev"
  )
  
  message(paste0("Fitting psoriasis NMA for PASI response with ", effects, 
                 " effects, ", cutpoints, " cutpoints, and ", 
                 if (baseline == "adjusted") "baseline adjustment" else 
                   "no baseline adjustment"))
  
  jags(
    data = data, parameters.to.save = params, inits = NULL, 
    model.file = filename, n.chains = 2, n.iter = 2000, n.burnin = 1000, n.thin = 1
  )
}




# models <- list(
#   fe_fez_u = pso_jags(pasi_jags, filename = "JAGS/fe_fez_u.jags", effects = "fixed", cutpoints = "fixed", baseline = "unadjusted"),
#   re_fez_u = pso_jags(pasi_jags, filename = "JAGS/re_fez_u.jags", effects = "random", cutpoints = "fixed", baseline = "unadjusted"),
#   fe_rez_u = pso_jags(pasi_jags, filename = "JAGS/fe_rez_u.jags", effects = "fixed", cutpoints = "random", baseline = "unadjusted"),
#   re_rez_u = pso_jags(pasi_jags, filename = "JAGS/re_rez_u.jags", effects = "random", cutpoints = "random", baseline = "unadjusted"),
#   fe_fez_a = pso_jags(pasi_jags, filename = "JAGS/fe_fez_a.jags", effects = "fixed", cutpoints = "fixed", baseline = "adjusted"),
#   re_fez_a = pso_jags(pasi_jags, filename = "JAGS/re_fez_a.jags", effects = "random", cutpoints = "fixed", baseline = "adjusted"),
#   fe_rez_a = pso_jags(pasi_jags, filename = "JAGS/fe_rez_a.jags", effects = "fixed", cutpoints = "random", baseline = "adjusted"),
#   re_rez_a = pso_jags(pasi_jags, filename = "JAGS/re_rez_a.jags", effects = "random", cutpoints = "random", baseline = "adjusted")
# )
# 
# nma_results(models$fe_fez_u, effects = "fixed", method = "standard") |> View()
# process_jags(models$fe_fez_a)
# 
# models |> 
#   lapply(\(x) {process_jags(x)$summary}) |> 
#   bind_rows(.id = "id") |> 
#   filter(!is.na(drug), ! drug %in% c("Phototherapy", "Mirikizumab")) |> 
#   mutate(drug = forcats::fct_reorder(drug, mean, .fun = base::mean)) |>
#   ggplot(aes(x = mean, y = drug, colour = id)) +
#   geom_pointrange(aes(xmin = `2.5%`, xmax = `97.5%`), 
#                   position = position_dodge(width = 0.7), shape = 15, size = 0.1) +
#   scale_colour_viridis_d(option = "turbo") +
#   theme_minimal() +
#   theme(legend.position = "top")
# ggsave("output/forest.png", height = 20, width = 7)
# 
# summaries <- lapply(models, process_jags)
# 
# lapply(summaries, \(x) {x$DIC})
# 
# traceplot(models$fe_fez_u)
