source("R/utils/jags-process.R")
source("R/utils/load-chars.R") 
source("R/utils/load-pasi-reg.R") # Load PASI data only for trials with chosen characteristics
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")

input <- c(pasi_jags, char_list)
which(is.na(input$r[,,1]) != is.na(input$bio), arr.ind = TRUE)
input$timepoint <- input$timepoint / 4

j <- list()

j$unadjusted <- pso_jags(input, niter = 10000, effects = "random", cutpoints = "trial")
j$timepoint <- pso_jags(input, niter = 10000, effects = "random", cutpoints = "trial", reg = c("timepoint"), reg_mean = c(4))
j$covar <- pso_jags(input, niter = 20000, effects = "random", cutpoints = "trial", reg = c("timepoint", "weight", "bio"), reg_mean = c(4, 9, 0.3))
j$baseline <- pso_jags(input, niter = 20000, effects = "random", cutpoints = "trial", reg = c("timepoint", "mu"), reg_mean = c(4, 0.7))
j$cov_baseline <- pso_jags(input, niter = 20000, effects = "random", cutpoints = "trial", reg = c("timepoint", "mu", "weight", "bio"), reg_mean = c(4, 0.7, 9, 0.3))

check_fits(j)

forests(j)

jp <- lapply(j, process_jags)

saveRDS(jp, "R/saved-results/metareg-rezi.rds")
