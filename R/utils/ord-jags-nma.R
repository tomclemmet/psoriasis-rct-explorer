# Multinomial ordered NMA code
library(R2jags)
source("R/utils/jags-process.R")

# Function to write and save JAGS code according to chosen arguments and run
# the analysis
pso_jags <- function(
    data, 
    filename = "JAGS/temp.jags",
    niter = 2000,
    effects = c("fixed", "random"),
    cutpoints = c("fixed", "treatment", "trial",  "arm"),
    class = c("independent", "exchangeable"),
    consistency = c("consistency", "ume"),
    reg = c(),
    reg_mean = c(),
    reg_scale = c(),
    reg_effects = c(), # common, exchangeable, or independent
    reg_level = c(), # treatment or class (only relevant for for exchangeable or independent models)
    output = c("processed", "jags")
) {
  effects <- match.arg(effects)
  cutpoints <- match.arg(cutpoints)
  class <- match.arg(class)
  consistency <- match.arg(consistency)
  output <- match.arg(output)

  # Remove redundant data based on arguments, avoiding JAGS warning
  if (class == "independent" & !any(reg_level == "class")) {
    data$cl <- NULL
    data$ncl <- NULL
  }
  
  # Set default regression effects model
  if (length(reg) > 0 && length(reg_effects) == 0) {
    reg_effects <- rep("common", length(reg))
  }
  if (length(reg) > 0 && length(reg_level) == 0) {
    reg_level <- rep("treatment", length(reg))
  }
  
  # Generic code setting up the conditional binomial likelihood
  setup <- "
model {
  # *** PROGRAM STARTS
  for(i in 1:ns){                                                               # LOOP THROUGH STUDIES
    w[i, 1] <- 0                                                                # adjustment for multi-arm trials is zero for control arm
    delta[i, 1] <- 0                                                            # treatment effect is zero for control arm
    mu[i] ~ dnorm(0, .01)                                                        # vague priors for all trial baselines
    for (k in 1:na[i]) {                                                        # LOOP THROUGH ARMS
      p[i, k, 1] <- 1                                                           # Pr(PASI >0)
      for (j in 1:(nc[i] - 1)) {                                                # LOOP THROUGH CATEGORIES
        r[i, k, j] ~ dbin(q[i, k, j], n[i, k, j])                               # binomial likelihood
        q.raw[i, k, j] <- 1 - (p[i, k, C[i, j + 1]] /                           # conditional probabilities
        max(p[i, k, C[i, j]], 1e-14)) 
        q[i, k, j] <- max(1e-12, min(1 - 1e-12, q.raw[i, k, j]))
        "
  
  # Linear predictor theta with trial-specific baseline and treatment effect
  theta <- "theta[i, k, j] <- mu[i] + delta[i, k] + "
  
  # Term for cutpoints depends on chosen model
  z <- "z[C[i, j + 1] - 1] "
  z_tx <- "zeta[t[i, k], C[i, j + 1] - 1] "
  z_trial <- "zeta[i, C[i, j + 1] - 1] "
  z_arm <- "zeta[i, k, C[i, j + 1] - 1] "
  
  # Terms for meta-regression data and coefficients
  if (length(reg) > 0) {
    metareg <- c()
    for (i in 1:length(reg)) {
      t_index <- if (reg_level[i] == "treatment") "[t[i, k]] " else "[cl[t[i, k]]]"
      c_index <- if (reg_level[i] == "treatment") "[t[i, 1]]" else "[cl[t[i, 1]]]"
      metareg[i] <- paste0(
        "+ (beta_", reg[i], t_index, " - beta_", reg[i], c_index, ")",
        " * (", 
        reg[i], if(reg[i] %in% c("mu", "timepoint")) "[i]" else "[i, k]", " - ", reg_mean[i] / reg_scale[i], ")"
      )
    }
    metareg <- paste0(metareg, collapse = " ")
  }
  
  # Generic deviance calculations
  deviance <- "
        rhat[i, k, j] <- q[i, k, j] * n[i, k, j]                                # predicted number events
        dv[i, k, j] <- 2 * (                                                    # Deviance contribution of each category  
          r[i, k, j] * 
            (log(max(r[i, k, j], 1e-10)) - log(max(rhat[i, k, j], 1e-10))) +            
            (n[i, k, j] - r[i, k, j]) * 
            (log(max(n[i, k, j] - r[i, k, j], 1e-10)) - 
               log(max(n[i, k, j] - rhat[i, k, j], 1e-10)))
        )
      }
      dev[i, k] <- sum(dv[i, k, 1:(nc[i] - 1)])                                 # deviance contribution of each arm"
  
  # Generic code linking theta to predicted probability of response
  phi <- "
      for (j in 2:nc[i]) {                                                      # LOOP THROUGH CATEGORIES
        p[i, k, C[i, j]] <- phi.adj[i, k, j]                                    # link function
        # adjust link function phi(x) for extreme values that can give numerical errors
        # when x< -8, phi(x)=0, when x> 8, phi(x)=1
        phi.adj[i, k, j] <- step(8 + theta[i, k, j - 1]) * 
          (step(theta[i, k, j - 1] - 8) +
             step(8 - theta[i, k, j - 1]) * phi(theta[i, k, j - 1]))
      }
    }"
  
  # Trial-specific delta terms for random effects model
  delta_re <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] ~ dnorm(md[i, k], taud[i, k])
      md[i, k] <- d[t[i, k]] - d[t[i, 1]] + sw[i, k]                            # mean of LHR distributions, with multi-arm trial correction
      taud[i, k] <- tau * 2 * (k - 1) / k                                       # precision of LHR distributions (with multi-arm trial correction)
      w[i, k] <- (delta[i, k] - d[t[i, k]] + d[t[i, 1]])                        # adjustment, multi-arm RCTs
      sw[i, k] <- sum(w[i, 1:(k-1)]) / (k-1)                                    # cumulative adjustment for multi-arm trials
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  # Trial-specific delta terms with no indirect evidence for UME model
  delta_re_ume <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] ~ dnorm(d[t[i, 1], t[i, k]], tau)
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  # Fixed delta terms for fixed effects model
  delta_fe <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] <- d[t[i, k]] - d[t[i, 1]]
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  # Fixed delta terms with no indirect evidence for UME model
  delta_fe_ume <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] <- d[t[i, 1], t[i, k]]
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  # Priors for different types of cutpoints
  fez <- "
  z[1] <- 0                                                                     # set z50=0
  for (j in 2:(Cmax-1)) {                                                       # Set priors for z, for any number of categories
    z.aux[j] ~ dunif(0,5)                                                       # priors
    z[j] <- z[j-1] + z.aux[j]                                                   # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
  }"
  rez_tx <- "
  for (i in 1:nt) {zeta[i, 1] <- 0 } # set z50=0
  z[1] <- 0
  for (j in 2:(Cmax-1)) {                                                       # Set priors for z, for any number of categories
    z.aux[j] ~ dunif(0,5)                                                       # priors
    z[j] <- z[j - 1] + z.aux[j]
    for (i in 1:nt) {
      zeta.aux[i, j] ~ dnorm(z.aux[j], tau_z)
      zeta[i, j] <- zeta[i, j - 1] + zeta.aux[i, j]                             # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
    }
  }"
  rez_trial <- "                                                                    
  for (i in 1:ns) { zeta[i, 1] <- 0 } # set z50=0
  z[1] <- 0
  for (j in 2:(Cmax-1)) {                                                       # Set priors for z, for any number of categories
    z.aux[j] ~ dunif(0,5)                                                       # priors
    z[j] <- z[j - 1] + z.aux[j]
    for (i in 1:ns) {
      zeta.aux[i, j] ~ dnorm(z.aux[j], tau_z)
      zeta[i, j] <- zeta[i, j - 1] + zeta.aux[i, j]                             # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
    }
  }"
  rez_arm <- "                                                                    
  for (i in 1:ns) {
    for (k in 1:na[i]) {
      zeta[i, k, 1] <- 0 # set z50=0
    }
  } 
  z[1] <- 0
  for (j in 2:(Cmax-1)) {                                                       # Set priors for z, for any number of categories
    z.aux[j] ~ dunif(0,5)                                                       # priors
    z[j] <- z[j - 1] + z.aux[j]
    for (i in 1:ns) {
      for (k in 1:na[i]) {
        zeta.aux[i, k, j] ~ dnorm(z.aux[j], tau_z)
        zeta[i, k, j] <- zeta[i, k, j - 1] + zeta.aux[i, k, j]                  # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
      }
    }
  }"
  
  # Treatment effect priors
  d_priors <- "
  d[1] <- 0                                                                     # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) }                                      # vague priors for treatment effects
"
  
  # Treatment effect priors for all comparisons in UME model
  d_priors_ume <- "
  for (k in 1:nt) { d[k,k] <- 0 }                                               # treatment effect is zero for reference treatment
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      d[c, k] ~ dnorm(0, 0.0001)
      d[k, c] <- -d[c, k]                                                       # allows treatments to be out of order
    }
  }
"
  # Class- and treatment-effect priors
  d_priors_class <- "
  d[1] <- 0
  m[1] <- 0
  for (k in 2:nt) { d[k] ~ dnorm(m[cl[k]], tau_cl) }
  for (p in 2:ncl){ m[p] ~ dnorm(0,.0001) }  
"
  
  if (length(reg) > 0) {
    reg_priors <- c()
    for (i in 1:length(reg)) {
      max <- if (reg_level[i] == "treatment") "nt" else "ncl"
      reg_priors[i] <- paste(
        paste0("\n  beta_", reg[i], "[1] <- 0"),
        if (reg_effects[i] == "common") {
          paste0("for (k in 2:", max, ") { beta_", reg[i], "[k] <- B_", reg[i], "}")
        } else if (reg_effects[i] == "exchangeable") {
          paste0("for (k in 2:", max, ") { beta_", reg[i], "[k] ~ dnorm(B_", reg[i], ", tau_", reg[i] ,")}
  sd_", reg[i], " ~ dunif(0, 5)
  tau_", reg[i], " <- pow(sd_", reg[i], ", -2)")
        } else if (reg_effects[i] == "independent") {
          paste0("for (k in 2:", max, ") { beta_", reg[i], "[k] ~ dnorm(0, .01) }")
        },
        paste0("B_", reg[i], "  ~ dnorm(0,.01)"),
        sep = "\n  "
      )
    }
    reg_priors <- paste0(reg_priors, collapse = "")
  }
  
  # Other priors (some may not be used, depending on model)
  priors <- "
  totresdev <- sum(resdev[1:ns])                                                    # Total Residual Deviance

  sd ~ dunif(0, 5)                                                              # vague prior for between-trial SD
  sd_z ~ dunif(0, 5)
  sd_cl ~ dunif(0, 5)
  tau <- pow(sd,-2)   
  tau_z <- pow(sd_z, -2) # between-trial precision = (1 / between-trial variance)
  tau_cl <- pow(sd_cl, -2)
  mubar <- mean(mu[])
  
  
  # A ~ dnorm(1.097,123) 
  p0 ~ dbeta(113.2, 566.2)
  A <- probit(1 - p0)
  # calculate prob of achieving PASI 50/75/90/100 on treatment k"

  # Probability calculations for fixed or random trial/arm-level cutpoints
  probs_fez <- "
  for (k in 1:nt) {
    for (j in 1:(Cmax - 1)) { prob[j,k] <- phi(A + d[k] + z[j]) }
  }"
  
  # Probability calculations for random treatment-level cutpoints
  probs_rez_tx <- "
  for (k in 1:nt) {
    for (j in 1:(Cmax - 1)) { prob[j,k] <- phi(A + d[k] + zeta[k, j]) }
  }"
  
  # End of file
  end <- "
  # *** PROGRAM ENDS 
}"
  
  # Combine relevant code segments into a single string
  model_code <- paste0(
    setup,
    theta,
    switch(cutpoints,
           "fixed" = z,
           "treatment" = z_tx,
           "trial" = z_trial,
           "arm" = z_arm),
    if (length(reg) > 0) metareg,
    deviance,
    phi,
    switch(paste(effects, consistency),
           "random consistency" = delta_re,
           "random ume" = delta_re_ume,
           "fixed consistency" = delta_fe,
           "fixed ume" = delta_fe_ume),
    switch(cutpoints,
           "fixed" = fez,
           "treatment" = rez_tx,
           "trial" = rez_trial,
           "arm" = rez_arm),
    if (consistency == "consistency") {if (class == "exchangeable") d_priors_class else d_priors} else d_priors_ume,
    if (length(reg) > 0) reg_priors,
    priors,
    if (consistency == "consistency") {if (cutpoints != "treatment") probs_fez else probs_rez_tx} else NULL,
    end
  )
  
  # # Choose filename based on arguments
  # if(is.na(filename)) {
  #   filename <- paste0("JAGS/", paste(
  #     if_else(effects == "random", "re", "fe"),
  #     switch(cutpoints,
  #            "fixed" = "fez",
  #            "treatment" = "rezt",
  #            "trial" = "rezi",
  #            "arm" = "reza"),
  #     if_else(length(reg) > 0, "a", "u"),
  #     if_else(class == "exchangeable", "c", "nc"),
  #     if_else(consistency == "ume", "ume", "con"),
  #     sep = "_"
  #   ), ".jags")
  # }
  
  # Write model code to file
  writeLines(model_code, filename)

  # Choose parameters of interest for JAGS to track
  params <- c(
    "d", "z", "mu",
    if (effects == "random") c("sd", "delta", "zeta") else NULL,
    if (cutpoints == "fixed") NULL else "sd_z",
    if (length(reg) > 0 & any(reg_effects != "independent")) c(paste0("B_", reg), "mubar") else NULL,
    if (any(reg_effects != "common")) paste0("beta_", reg) else NULL,
    if (class == "exchangeable") c("m", "sd_cl") else NULL,
    if (consistency == "ume") NULL else "prob",
    "totresdev",
    "dv", "rhat"
  )

  # Specify reasonable initial values for each chain
  inits <- list(
    list(
      d = c(NA, rep(0, data$nt - 1)),
      mu = rep(0, data$ns),
      z.aux = c(NA, rep(1, data$Cmax-2))
    ),
    list(
      d = c(NA, rep(1, data$nt - 1)),
      mu = rep(0, data$ns),
      z.aux = c(NA, rep(1, data$Cmax-2))
    )
  )
  if (effects == "random") {
    inits[[1]]$sd <- 1
    inits[[2]]$sd <- 0.5
  }
  if (cutpoints != "fixed") {
    inits[[1]]$sd_z <- 1
    inits[[2]]$sd_z <- 0.5
  }
  if (class == "exchangeable") {
    inits[[1]]$m <- c(NA, rep(0, data$ncl - 1))
    inits[[1]]$sd_cl <- 1
    inits[[2]]$m <- c(NA, rep(1, data$ncl - 1))
    inits[[2]]$sd_cl <- 0.5
  }
  if (consistency == "ume") {
    inits[[1]]$d <- NULL
    inits[[2]]$d <- NULL
  }
  if (length(reg) > 0) {
    for (i in 1:length(reg)) {
      inits[[1]][[paste0("B_", reg[i])]] <- 0
      inits[[2]][[paste0("B_", reg[i])]] <- 0
    }
  }

  # Scale and centre parameters
  if (length(reg) > 0) {
    for (i in 1:length(reg)) {
      data[[reg[i]]] <- data[[reg[i]]] / reg_scale[i]
    }
  }
  
  # Send message with current model type
  message(paste0("Fitting ", if_else(consistency == "ume", "UME ", ""), "NMA for PASI response with ", effects,
                 " effects, ", cutpoints, if_else(cutpoints == "fixed", "", "-level"), " cutpoints, ",
                 if (length(reg) > 0) paste0("meta-regression on ", paste(reg, collapse = ", ")),
                 if (class == "exchangeable") ", exchangeable class effects" else NULL))

  # Run model
  set.seed(123)
  fit <- jags(
    data = data, parameters.to.save = params, inits = inits,
    model.file = filename, n.chains = 2, n.iter = niter, n.burnin = niter/2, n.thin = ceiling(niter/20000)
  )

  # Send informative warning if any parameters have failed to converge
  Rhat <- fit$BUGSoutput$summary[, "Rhat"]
  if (any(Rhat > 1.1)) {
    warning(paste0(
      "WARNING: The following parameters have an Rhat greater than 1.1: ",
      paste(paste0(names(Rhat)[Rhat > 1.1], "(", round(Rhat[Rhat > 1.1], 3), ")"), collapse = ", ")
    ))
  }

  fit
}
