# Bivariate NMA code based on approach in TSD20 appendix D.2 (p.116)
rm(list = ls())
library(DBI)
library(dplyr)
library(tidyr)
library(R2jags)


drug_order = c(
    "Placebo", "Acitretin", "Adalimumab", "Apremilast", "Bimekizumab",
    "Brodalumab", "Certolizumab", "Cyclosporin", "Deucravacitinib", "Etanercept",
    "Fumaric acid esters", "Guselkumab", "Icotrokinra", "Infliximab", "Ixekizumab",
    "Izokibep", "Methotrexate", "Mirikizumab", "Netakimab", "Orismilast",
    "Phototherapy", "Risankizumab", "Roflumilast", "Secukinumab", "Sonelokimab",
    "Tildrakizumab", "Tofacitinib", "Ustekinumab", "Xeligekimab", "Zasocitinib"
  )

con <- dbConnect(RSQLite::SQLite(), "app/psoriasis-rcts.sqlite")

pasi <- dbReadTable(con, "v_pasi") |> 
  filter(!is.na(pasi75)) |> 
  select(trial, ref_id, drug, n, pasi75) |> 
  summarise(.by = c(ref_id, drug), n_pasi = sum(n), pasi75 = sum(pasi75)) |> 
  mutate(t_pasi = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))))


dlqi <- dbReadTable(con, "v_dlqi") |> 
  filter(!is.na(dlqi_0_1)) |> 
  select(trial, ref_id, drug, n, dlqi_0_1) |> 
  summarise(.by = c(ref_id, drug), n_dlqi = sum(n), dlqi_0_1 = sum(dlqi_0_1)) |> 
  mutate(t_dlqi = as.numeric(factor(drug, levels = c("Placebo", setdiff(sort(drug), "Placebo")))))
dbDisconnect(con)

t_index <- data.frame(drug = drug_order) |> 
  left_join(distinct(pasi, drug, t_pasi), by = "drug") |> 
  left_join(distinct(dlqi, drug, t_dlqi), by = "drug")

bivar_data <- full_join(pasi, dlqi, by = c("ref_id", "drug")) |> 
  mutate(.by = ref_id, arm_no = row_number(), na = n(),
         o1 = as.numeric(!is.na(n_pasi)), o2 = as.numeric(!is.na(n_pasi)) * 2) |> 
  select(-t_pasi, -t_dlqi) |> 
  left_join(t_index, by = "drug")

df2 <- bivar_data |> 
  select(ref_id, na, arm_no, t_pasi, t_dlqi, o1, o2) |> 
  pivot_longer(cols = starts_with("t_"), names_to = "outcome",
               names_prefix = "t_", values_to = "index") |> 
  pivot_wider(names_from = arm_no, names_prefix = "t", values_from = index)

df3 <- bivar_data |> 
  select(-c(drug, na:o2)) |> 
  mutate(y1 = pasi75 / n_pasi, y2 = dlqi_0_1 / n_dlqi, 
         se1 = sqrt((y1 / (1 - y1) / n_pasi)), se2 = sqrt((y2 / (1 - y2) / n_dlqi)), 
         v = NA, .keep = "unused") |> 
  arrange(!is.na(y2))

bivar_jags <- list(
  ns = n_distinct(bivar_data$ref_id),
  no = 2,
  nobs1 = sum(is.na(bivar_data$pasi75) | is.na(bivar_data$dlqi_0_1)),
  nobs2 = sum(!is.na(bivar_data$pasi75) & !is.na(bivar_data$dlqi_0_1)),
  nt.total = c(n_distinct(pasi$drug), n_distinct(dlqi$drug)),
  na1 = summarise(bivar_data, .by = ref_id, n = n())$n,
  
  s = df2$ref_id,
  t = select(df2, t1:t3) |> as.matrix(),
  na2 = df2$na,
  o = select(df2, o1:o2) |> as.matrix(),
  
  study = df3$ref_id,
  arm = df3$arm_no,
  y = select(df3, y1:y2) |> as.matrix(),
  se = select(df3, se1:se2) |> as.matrix(),
  v = df3$v
)

bivar_code <- function() {

  # Likelihood for studies reporting a single outcome (arm level data)
  for(i in 1:nobs1){
    tmp[i] <- v[i] #dummy variable within-study corr in datafile3
    omega1[i,1] <- pow(se[i,1],-2) #precision
    y[i,1] ~ dnorm(mean.y[study[i],arm[i],1],omega1[i,1]) # Normal dist.
  }
  # Likelihood for studies reporting two outcomes (arm level data )
  rhoW ~ dunif(-1,1) # uniform(a,b) prior distr. for within-study corr.
  for(i in (nobs1+1):(nobs1+nobs2)){
    omega2[i,1:no,1:no] <- inverse(cov.mat[i,,])
    y[i,1:no] ~ dmnorm(mean.y[study[i],arm[i],1:no],omega2[i,,]) # mvNorm distr.
    # Define within-study covariance matrix
    cov.mat[i,1,1] <- pow(se[i,1],2)
    cov.mat[i,2,2] <- pow(se[i,2],2)
    cov.mat[i,1,2] <- se[i,1]*se[i,2] *rhoW
    cov.mat[i,2,1] <- cov.mat[i,1,2]
  }
  # Transform unobserved effects to one-row per study and take contrasts
  for(j in 1:ns) {
    for(k in 1:na1[j]) {
      for(m in 1:no) {
        mean.y[j,k,m] <- mu[j,m] + delta[j,k,m]
      }
    }
  }
  # Multivariate random-effects (homogenous variance model)
  for(j in 1:ns){
    for(m in 1:no) {
      delta[j,1,m] <-0 #set delta in treatment 1/baseline to zero
      w[j,1,m] <-0 #set multi-arm adjustment in trt 1 to zero
    }
    for(k in 2:na1[j]){
      delta[j,k,1:no] ~ dmnorm(md[j,k,1:no],PREC[j,k,1:no,1:no])
      #Set precision matrix T for multiple outcomes
      for(m in 1:no) {
        for(mm in 1:no) {
          PREC[j,k,m,mm] <- T[m,mm]*2*(k-1)/k
        }
      }
    }
  }
  # Consistency relations between basic parameters
  for(i in 1:(ns*no)) {
    for(k in 2:na2[i]) {
      md[s[i],k,o[i,1]] <- (d[o[i,1],t[i,k]] - d[o[i,1],t[i,1]]) * 
                            equals(o[i,1],o[i,2]) + sw[s[i],k,o[i,1]]
      w[s[i],k,o[i,1]] <- (delta[s[i],k,o[i,1]] - 
                          (d[o[i,1],t[i,k]] - d[o[i,1],t[i,1]])) * 
                          equals(o[i,1],o[i,2])
      sw[s[i],k,o[i,1]] <- sum(w[s[i],1:k-1,o[i,1]]) / (k-1)
    }
  }
  ##########################################################
  # Constraints.
  # set effect in trt 1 on both outcome to zero.
  # outcome 2 has 8 treatments buts outcome 2 has only 6 treatments#
  # so need to set d[2,7] and d[2,1]
  d[1,1] <-0
  d[2,1] <-0
  d[2,7] <-0
  d[2,8] <-0

  # trt effects exponentiated and prior distributions
  for(m in 1:no) {
    sd[m] ~ dunif(0, 5)
    sigma[m,m] <- pow(sd[m],2)
    for(k in 2:nt.total[m]){
      or[m,k] <- exp(d[m,k])
      d[m,k] ~ dnorm(0,0.0001)
    }
  }
  # Prior distributions and estimated between study correlation matrix
  for(j in 1:ns){
    for(m in 1:no) {
      mu[j, m] ~ dnorm(0,0.0001)
    }
  }
  # Parameterization of the between-studies covariance based
  # on the separation strategy by spherical decomposition (Wei et al 2013)
  T[1:no,1:no] <- inverse(sigma[,])
  pi <- 3.1415
  for(i in 1:2) {
    for(j in (i+1):no) {
      sigma[i,j] <- rhoB[i,j]*sd[i]*sd[j]
      sigma[j,i] <- sigma[i,j]
      g[j,i] <- 0
      a[i,j] ~ dunif(0, pi)
      rhoB[i,j] <- inprod(g[,i], g[,j])
    }
  }
  g[1,1] <- 1
  g[1,2] <- cos(a[1,2])
  g[2,2] <- sin(a[1,2])
  # Additional parameterization in 3-outcome problem
  # g[1,3] <- cos(a[1,3])
  # g[2,3] <- sin(a[1,3])*cos(a[2,3])
  # g[3,3] <- sin(a[1,3])*sin(a[2,3])
  # pairwise ORs and LORs for all possible pair-wise comparisons, if nt>2
  # outcome 1
  for (c in 1:(nt.total[1]-1)) {
    for (k in (c+1):nt.total[1]) {
      or2[c,k] <- exp(d[1,k] - d[1,c])
      lor1[c,k] <- (d[1,k]-d[1,c])
    }
  }
  # outcome 2
  for (c in 1:(nt.total[2]-1)) {
    for (k in (c+1):nt.total[2]) {
      or2[c,k] <- exp(d[2,k] - d[2,c])
      lor2[c,k] <- (d[2,k]-d[2,c])
    }
  }
}
lapply(bivar_jags, head)
m <- jags(
  data = bivar_jags, parameters.to.save = c("d"),
  inits = NULL, model.file = bivar_code, n.chains = 2, n.iter = 1000,
  n.burnin = 500, n.thin = 1
)
