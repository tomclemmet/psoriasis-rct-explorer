library(R2jags)

pso_jags <- function(
    data, 
    filename = NA,
    niter = 2000,
    effects = c("fixed", "random"),
    cutpoints = c("fixed", "treatment", "trial",  "arm"),
    baseline = c("unadjusted", "adjusted"),
    class = c("independent", "exchangeable"),
    consistency = c("consistency", "ume")
) {
  effects <- match.arg(effects)
  cutpoints <- match.arg(cutpoints)
  baseline <- match.arg(baseline)
  class <- match.arg(class)
  consistency <- match.arg(consistency)
  
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
        max(p[i, k, C[i, j]], 1e-14)) 
        q[i, k, j] <- max(1e-12, min(1 - 1e-12, q.raw[i, k, j]))
        "
  theta <- "theta[i, k, j] <- mu[i] - delta[i, k] + "
  
  z <- "z[C[i, j + 1] - 1]"
  z_tx <- "zeta[t[i, k], C[i, j + 1] - 1]"
  z_trial <- "zeta[i, C[i, j + 1] - 1]"
  z_arm <- "zeta[i, k, C[i, j + 1] - 1]"
  
  baseline_adj <- " + beta * (mu[i] - mmu) * (1 - equals(k, 1))"
  
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
  
  phi <- "
      for (j in 2:nc[i]) {                                                      # LOOP THROUGH CATEGORIES
        p[i, k, C[i, j]] <- 1 - phi.adj[i, k, j]                                # link function
        # adjust link function phi(x) for extreme values that can give numerical errors
        # when x< -8, phi(x)=0, when x> 8, phi(x)=1
        phi.adj[i, k, j] <- step(8 + theta[i, k, j - 1]) * 
          (step(theta[i, k, j - 1] - 8) +
             step(8 - theta[i, k, j - 1]) * phi(theta[i, k, j - 1]))
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
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  delta_re_ume <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] ~ dnorm(d[t[i, 1], t[i, k]], tau)
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  delta_fe <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] <- d[t[i, k]] - d[t[i, 1]]
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
  delta_fe_ume <- "
    for (k in 2:na[i]) {                                                        # LOOP THROUGH ARMS
      delta[i, k] <- d[t[i, 1], t[i, k]]
    }
    resdev[i] <- sum(dev[i, 1:na[i]])                                           # summed residual deviance contribution for this trial
  }"
  
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
      zeta.aux[i, j] ~ dnorm(z.aux[j], tauz)
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
      zeta.aux[i, j] ~ dnorm(z.aux[j], tauz)
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
        zeta.aux[i, k, j] ~ dnorm(z.aux[j], tauz)
        zeta[i, k, j] <- zeta[i, k, j - 1] + zeta.aux[i, k, j]                  # ensures z[j]~Uniform(z[j-1], z[j-1]+5)
      }
    }
  }"
  
  d_priors <- "
  d[1] <- 0                                                                     # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) }                                      # vague priors for treatment effects
"
  d_priors_ume <- "
  for (k in 1:nt) { d[k,k] <- 0 }                                               # treatment effect is zero for reference treatment
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      d[c, k] ~ dnorm(0, 0.0001)
      d[k, c] <- -d[c, k]                                                       # allows treatments to be out of order
    }
  }
"
  d_priors_class <- "
  d[1] <- 0
  m[1] <- 0
  for (k in 2:nt) {
    d[k] ~ dnorm(m[cl[k]], taucl)
  }
  for (p in 2:ncl){ m[p] ~ dnorm(0,.0001) }  
"
  
  priors <- "
  totresdev <- sum(resdev[])                                                    # Total Residual Deviance
  beta ~ dnorm(0,.0001)
  sd ~ dunif(0, 5)                                                              # vague prior for between-trial SD
  sdz ~ dunif(0, 5)
  sdcl ~ dunif(0, 5)
  tau <- pow(sd,-2)   
  tauz <- pow(sdz, -2) # between-trial precision = (1 / between-trial variance)
  taucl <- pow(sdcl, -2)
  mubar <- mean(mu[])
  
  
  # A ~ dnorm(1.097,123) 
  p0 ~ dbeta(113.2, 566.2)
  A <- probit(1 - p0)
  # calculate prob of achieving PASI 50/75/90/100 on treatment k"

  probs_fez <- "
  for (k in 1:nt) {
    for (j in 1:(Cmax - 1)) { prob[j,k] <- 1 - phi(A - d[k] + z[j]) }
  }"
  
  probs_rez_tx <- "
  for (k in 1:nt) {
    for (j in 1:(Cmax - 1)) { prob[j,k] <- 1 - phi(A - d[k] + zeta[k, j]) }
  }"
  
  end <- "
  # *** PROGRAM ENDS 
}"
  
  model_code <- paste0(
    setup,
    theta,
    switch(cutpoints,
           "fixed" = z,
           "treatment" = z_tx,
           "trial" = z_trial,
           "arm" = z_arm),
    if (baseline == "unadjusted") "" else baseline_adj,
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
    priors,
    if (consistency == "consistency") {if (cutpoints != "treatment") probs_fez else probs_rez_tx} else NULL,
    end
  )
  
  if(is.na(filename)) {
    filename <- paste0("JAGS/", paste(
      if_else(effects == "random", "re", "fe"),
      switch(cutpoints,
             "fixed" = "fez",
             "treatment" = "rezt",
             "trial" = "rezi",
             "arm" = "reza"),
      if_else(baseline == "adjusted", "a", "u"),
      if_else(class == "exchangeable", "c", "nc"),
      if_else(consistency == "ume", "ume", "con"),
      sep = "_"
    ), ".jags")
  }
  
  writeLines(model_code, filename)
  
  params <- c(
    "d", "z", "prob",
    if (effects == "random") "sd" else NULL,
    if (cutpoints == "fixed") NULL else "sdz",
    if (baseline == "adjusted") c("beta", "mubar") else NULL,
    if (class == "exchangeable") c("m", "sdcl") else NULL,
    "totresdev",
    "dv", "rhat"
  )
  
  message(paste0("Fitting ", if_else(consistency == "ume", "UME ", ""), "NMA for PASI response with ", effects, 
                 " effects, ", cutpoints, if_else(cutpoints == "fixed", "", "-level"), " cutpoints, ", 
                 if_else(baseline == "adjusted", "baseline adjustment", "no baseline adjustment"),
                 if (class == "exchangeable") ", exchangeable class effects" else NULL))
  
  jags(
    data = data, parameters.to.save = params, inits = NULL, 
    model.file = filename, n.chains = 2, n.iter = niter, n.burnin = niter/2, n.thin = 1
  )
}
