rm(list = ls())
source("R/meta-analyse/ma-utils.R")
source("R/meta-analyse/wide_format.R")
source("R/meta-analyse/pasi-jags-nma.R")
j <- readRDS("R/meta-analyse/jags_fits.rds")

j$re_rez_u_nc_con <- pso_jags(pasi_jags, effects = "random", cutpoints = "random")
j$re_rez_u_nc_ume <- pso_jags(pasi_jags, effects = "random", cutpoints = "random",
                              consistency = "ume")

process_jags(j$re_rez_u_nc_con)$dev_table |> View()

devplot(j$re_rez_u_nc_con, j$re_rez_u_nc_ume)

devplot(j$re_rez_u_nc_con, j$re_rez_u_nc_ume, "table") |> View()

# saveRDS(j, "R/meta-analyse/jags_fits.rds")