t1 <- "Risankizumab"
t2 <- "Secukinumab"

direct_data <- gen_pasi_jags(c(t1, t2))
if (nrow(pasi_wide) == 1) {
  direct_data <- list(
    nc = pasi_wide$nc,
    C = direct_data$C |> as.vector(),
    r = rbind(
      select(pasi_wide, starts_with("a1r")) |> as.vector(),
      select(pasi_wide, starts_with("a2r")) |> as.vector()
    ),
    n = rbind(
      select(pasi_wide, starts_with("a1n")) |> as.vector(),
      select(pasi_wide, starts_with("a2n")) |> as.vector()
    )
  )
  direct <- jags(
    data = direct_data,
    parameters.to.save = "d", inits = NULL, model.file = "JAGS/single_study.jags", 
    n.chains = 2, n.iter = 100000, n.burnin = 50000
  )
  dir <- posterior::as_draws_df(direct$BUGSoutput$sims.array)$d
} else {
  direct <- pso_jags(direct_data, effects = "random", cutpoints = "trial", baseline = "adjusted")
  dir <- posterior::as_draws_df(direct$BUGSoutput$sims.array)$`d[2]`
}

indirect_data <- gen_pasi_jags(c(t1, t2), direct = FALSE)
indirect <- pso_jags(indirect_data, effects = "random", cutpoints = "trial")

p1 <- paste0("d[", match(t1, pasi_drugs), "]")
p2 <- paste0("d[", match(t2, pasi_drugs), "]")

indirect_trace <- posterior::as_draws_df(indirect$BUGSoutput$sims.array)

data.frame(t1 = indirect_trace[[p1]], t2 = indirect_trace[[p2]]) |> 
  mutate(indir = t2 - t1, dir = dir) |> 
  select(-t1, -t2) |> 
  pivot_longer(everything(), names_to = "lab", values_to = "d") |> 
  ggplot() +
  geom_density(aes(d, colour = lab, fill = lab), alpha = 0.5)