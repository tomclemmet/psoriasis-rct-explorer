input <- c(pasi_jags, char_list)
which(is.na(input$r[,,1]) != is.na(input$bio), arr.ind = TRUE)
input$timepoint <- input$timepoint / 4

m <- pso_jags(
  input, niter = 10000,
  reg = c("timepoint", "age", "sex", "weight", "bio"), 
  reg_mean = c(4, 4, 0.3, 0.5, 8, 1.5, 0.2),
  effects = "random",
  cutpoints = "trial"
)

m2 <- pso_jags(
  input, niter = 10000,
  reg = c("timepoint", "bio"), 
  reg_mean = c(1.6, 0.2),
  effects = "random",
  cutpoints = "trial"
)


process_jags(m2)
forest(m2)
traceplot(m2, varname = "B_")

input$bio |> View()

pasi_wide$ref_id |> paste(collapse = ", ")

pasi_condensed |> filter(ref_id %in% c)