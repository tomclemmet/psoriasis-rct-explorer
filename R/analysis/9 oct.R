options(digits = 3)
source("R/utils/jags-process.R")
source("R/utils/ord-jags-nma.R") 
source("R/utils/jags-process.R")
source("R/utils/plots.R")
source("R/utils/load-chars.R") 
source("R/utils/load-pasi.R")
theme_set(theme_classic())

j <- list()

pasi_jags$timepoint <- log(pasi_jags$timepoint, 2)
j$rezi_u <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial"
)

j$rezi_u_t <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "independent", reg_level = "class"
)

# pasi_long |> 
#   distinct(trial) |> 
#   mutate(yr = if_else(
#     str_sub(trial, -4, -1)) |> 
#   arrange(yr)
j$rezi_u_t_yr <- pso_jags(
  pasi_jags, niter = 10000, effects = "random", cutpoints = "trial",
  reg = "timepoint", reg_mean = 16, reg_scale = 4, reg_effects = "independent", reg_level = "class",
)

bio <- load_chars("bio")
j$rezi_u_t_bio <- pso_jags(
  bio, niter = 10000, effects = "random", cutpoints = "trial",
  reg = c("timepoint", "bio"), reg_mean = c(16, 0.25), reg_scale = c(4, 1), 
  reg_effects = c("independent", "independent"), reg_level = c("class", "treatment")
)
jc$rezi_u_t_bio <- process_jags(j$rezi_u_t_bio)
forests(process_jags(j$rezi_u_t_bio), grp = "Regression coefficient") +
  xlim(-15, 15)
# jc <- lapply(j, process_jags)

forests(jc)
check_fits(jc)

pasi_long |> 
  left_join(pasi_id_lookup, by = c("drug", "t")) |> 
  filter(class == "conventional") |> 
  mutate(p = pasi75 / n) |> 
  ggplot(aes(x = timepoint, y = p)) +
  geom_point(aes(size = n, colour = drug), alpha = 0.7) +
  geom_smooth(aes(weight = n, colour = drug), method = "lm", se = FALSE) +
  geom_smooth(aes(weight = n), method = "lm", se = FALSE, colour = "black") +
  scale_x_continuous(breaks = seq(12, 24, 4))

pasi_long |> 
  filter(drug=="Fumaric acid esters", timepoint==24) |> 
  select(-pop_res, -pasi_high_rob) |> 
  mutate(p = pasi75/n) |> arrange(p)
  

bio <- load_chars("bio")

pasi_long |> 
  left_join(pasi_id_lookup, by = c("drug", "t")) |> 
  mutate(bio = as.vector(t(bio$bio))[!is.na(as.vector(t(bio$bio)))]) |>
  filter(class == "il12_23") |> 
  mutate(p = pasi75 / n) |> 
  ggplot(aes(x = bio, y = p)) +
  geom_point(aes(size = n, colour = drug), alpha = 0.7) +
  geom_text(aes(label = ref_id)) +
  geom_smooth(aes(weight = n, colour = drug), method = "lm", se = FALSE)# +
  # geom_smooth(aes(weight = n), method = "lm", colour = "black")

pasi_long |> 
  left_join(pasi_id_lookup, by = c("drug", "t")) |> 
  mutate(bio = as.vector(t(bio$bio))[!is.na(as.vector(t(bio$bio)))]) |>
  mutate(p = pasi100 / n) |> select(-c(t:timepoint)) |> 
  filter(drug == "Ustekinumab") |> arrange(bio)
