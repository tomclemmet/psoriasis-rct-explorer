source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R") 

j <- list()
jp <- list()

jp <- readRDS("R/saved-results/metareg-rezi.rds") 

# Unadjusted model
source("R/utils/load-pasi.R")
j$unadjusted_rezi <- pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4
)
jp$unadjusted_rezi <- process_jags(j$unadjusted_rezi)

mu_rezi <- pso_jags(
  pasi_jags, niter = 20000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "mu"), reg_mean = c(16, 0.7), reg_scale = c(4, 1)
)
mu_rezi

forests(jp$unadjusted_rezi, mu_rezi)

age <- load_chars("age")
j$age <- pso_jags(
  age, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "age"), reg_mean = c(16, 45), reg_scale = c(4, 10)
)
jp$age <- process_jags(j$age)


sex <- load_chars("male")
j$sex <- pso_jags(
  sex, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "male"), reg_mean = c(16, 0.7), reg_scale = c(4, 1)
)
jp$sex <- process_jags(j$sex)


ethnicity <- load_chars("white")
j$ethnicity <- pso_jags(
  ethnicity, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "white"), reg_mean = c(16, 0.9), reg_scale = c(4, 1)
)
jp$ethnicity <- process_jags(j$ethnicity)

weight <- load_chars("weight")
j$weight <- pso_jags(
  weight, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "weight"), reg_mean = c(16, 90), reg_scale = c(4, 10)
)
jp$weight <- process_jags(j$weight)


bmi <- load_chars("bmi")
j$bmi <- pso_jags(
  bmi, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "bmi"), reg_mean = c(16, 30), reg_scale = c(4, 5)
)
jp$bmi <- process_jags(j$bmi)


pasi <- load_chars("pasi")
j$pasi <- pso_jags(
  pasi, niter = 20000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "pasi"), reg_mean = c(16, 20), reg_scale = c(4, 10)
)
jp$pasi <- process_jags(j$pasi)


dur <- load_chars("duration")
j$dur <- pso_jags(
  dur, niter = 20000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "duration"), reg_mean = c(16, 15), reg_scale = c(4, 10)
)
jp$dur <- process_jags(j$dur)

psa <- load_chars("psa")
j$psa <- pso_jags(
  psa, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "psa"), reg_mean = c(16, 0.2), reg_scale = c(4, 1)
)
jp$psa <- process_jags(j$psa)

bio <- load_chars("bio")
j$bio <- pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "bio"), reg_mean = c(16, 0.2), reg_scale = c(4, 1)
)
jp$bio <- process_jags(j$bio)

j$bio <- pso_jags(
  bio, niter = 20000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "bio"), reg_mean = c(16, 0.2), reg_scale = c(4, 1), reg_effects = c("common", "independent")
)
jp$bio2 <- process_jags(j$bio)
forest(jp$bio2)

source("R/utils/load-pasi.R")

m0 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial"
)

m1 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "common", reg_level = "class"
)

m2 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "exchangeable", reg_level = "class"
)

m3 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "independent", reg_level = "class"
)

m4 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "exchangeable", reg_level = "treatment"
)


timepoint_mods <- list(
  unadjusted = process_jags(m0),
  common = process_jags(m1),
  ex_class = process_jags(m2),
  ind_class = process_jags(m3)
)

timepoint_mods$ex_class <- process_jags(m2)
timepoint_mods$ind_class <- process_jags(m3)
check_fits(timepoint_mods)
forests(timepoint_mods, grp = c("Heterogeneity", "Regression coefficient"))


ume0 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial",
  consistency = "ume"
)

ume1 <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "common", reg_level = "class",
  consistency = "ume"
)

ume2 <- pso_jags(
  pasi_jags, niter = 20000, effects = "random", cutpoints = "trial", 
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "exchangeable", reg_level = "class",
  consistency = "ume"
)

devplot(m0, ume0)
devplot(m1, ume1)
devplot(m2, ume2)
check_fits(list(m2, ume2))

devplot(m1, m2, output="table") |> arrang

check_fits(list(m0, m1, m2))
forests(m0, m1, m2, m3)
forests(m0, m1, m2, m3, grp = "Regression coefficient")
forests(m0, m1, m2, m3, grp = "Heterogeneity")

devplot(m0, m2, output="table") |> arrange(-diff)

timepoint_mods <- list(
  unadjusted = process_jags(m0),
  common = process_jags(m1),
  ex_class = process_jags(m2),
  ind_class = process_jags(m3),
  ex_tx = process_jags(m4)
)
# timepoint_mods <- readRDS("R/saved-results/timepoint_mods.rds")
check_fits(timepoint_mods)
forests(m0, m1, m2, m3)
saveRDS(timepoint_mods, "R/saved-results/timepoint_mods.rds")

bio <- load_chars("bio")
m0 <- pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4), 
  reg_effects = c("exchangeable"), reg_level = c("class")
)
m <- pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "bio"), reg_mean = c(16, 0.25), reg_scale = c(4, 1), 
  reg_effects = c("exchangeable", "common"), reg_level = c("class", "")
)
process_jags(m)
process_jags(m0)

# look into bio data

data.frame(
  b = as.vector(bio$bio),
  r = rep(1:78, 4),
  a = rep(1:4, 78)
) |> 
  filter(!is.na(t)) |> 
  left_join(pasi_id_lookup) |> 
  left_join(select(pasi_long)



forests(jp, grp = "B") +
  labs(title = "Univariate meta-regression", x = "Mean (95% CI)", y = "") +
  theme(legend.position = "none")
ggsave("output/metareg/B.png", height = 6, width = 5)

forests(jp, grp = "sd") +
  labs(title = "Univariate meta-regression", x = "Mean (95% CI)", y = "")
ggsave("output/metareg/sd.png", height = 6, width = 6)

jp_bio <- list()
bio <- load_chars("bio")
jp_bio$unadjusted <- process_jags(pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4)
))
jp_bio$bio <- process_jags(pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "bio"), reg_mean = c(16, 0.2), reg_scale = c(4, 1)
))
jp_bio$mu <- process_jags(pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "mu"), reg_mean = c(16, 0.7), reg_scale = c(4, 1)
))

check_fits(jp_bio) # Bio adjustment slightly reduces totresdev but increased DIC by ~ 5
devplot(jp_bio$bio, jp_bio$unadjusted)

ume_bio <- process_jags(pso_jags(
  bio, niter = 20000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint", "bio"), reg_mean = c(16, 0.2), reg_scale = c(4, 1),
  consistency = "ume"
))
devplot(jp_bio$bio, ume_bio)

ume <- process_jags(pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial", 
  reg = c("timepoint"), reg_mean = c(16), reg_scale = c(4),
  consistency = "ume"
))
devplot(jp_bio$unadjusted, ume)

bio <- load_chars("bio")
