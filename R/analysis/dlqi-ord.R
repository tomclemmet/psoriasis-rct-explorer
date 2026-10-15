options(digits = 3)
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
# source("R/utils/load-dlqi.R")
source("R/utils/load-chars-dlqi.R")

# j <- list()
# j$fe <- pso_jags(dlqi_jags, niter = 5000)
# j$re <- pso_jags(dlqi_jags, niter = 10000, effects = "random")
# j$rerezt <- pso_jags(dlqi_jags, niter = 10000, effects = "random", cutpoints = "treatment")
# j$rerezi <- pso_jags(dlqi_jags, niter = 10000, effects = "random", cutpoints = "trial")
# j$fea <- pso_jags(dlqi_jags, niter = 5000, reg = "mu", reg_mean = 1.3, reg_scale = 1)

check_fits(j)
forests(j)
compare_jags(j)

process_jags(j$fea)


bio <- load_chars("bio")

m <- list()
m$bio <- pso_jags(bio, niter = 20000, reg = "bio", reg_mean = 0.25, reg_scale = 1)
process_jags(m$bio)




