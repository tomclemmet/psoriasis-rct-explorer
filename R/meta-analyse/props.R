library(ggplot2)
library(dplyr)
library(stringr)
source("R/meta-analyse/wide_format.R")


j <- readRDS("R/meta-analyse/jags_fits.rds")

props_plot <- function(m) {
  
  outcomes <- c(
    "pasi50", "pasi75", "pasi90", "pasi100"
  )
  
  drug_rank <- m$results |> 
    filter(str_starts(param, "prob")) |> 
    mutate(drug = pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))]) |> 
    slice_head(n = 1, by = drug) |> 
    arrange(mean)
  
  m$results |> 
    filter(str_starts(param, "prob")) |> 
    mutate(
      drug = factor(
        pasi_drugs[as.numeric(str_extract(param, "(?<=,).*?(?=])"))],
        levels = drug_rank$drug
      ),
      outcome = factor(
        outcomes[as.numeric(str_extract(param, "(?<=\\[).*?(?=,)"))],
        levels = c("pasi50", "pasi75", "pasi90", "pasi100"),
        labels = c("PASI 50-75", "PASI 75-90", "PASI 90-100", "PASI 100")
      )
    ) |> 
    arrange(drug, desc(outcome)) |> 
    mutate(.by = drug, mean = mean - lag(mean, default = 0), .after = mean) |> 
    mutate(outcome = forcats::fct_rev(outcome)) |> 
    filter(drug %notin% c(
      "Icotrokinra", "Mirikizumab", "Netakimab", "Orismilast", "Roflumilast", 
      "Phototherapy", "Sonelokimab", "Tofacitinib", "Xeligekimab", "Zasocitinib")
    ) |> 
    ggplot(aes(x = mean, y = drug)) +
    geom_col(aes(fill = outcome), position = position_stack(reverse = TRUE)) +
    theme_classic() +
    scale_fill_viridis_d() +
    theme(legend.position = "right", axis.title.y = element_blank(), legend.title = element_blank(),
          plot.title = element_text(face = "bold")) +
    labs(title = "Predicted PASI Response", x = "Proportion")
}

props_plot(j$re_rezt_a_nc)
props_plot(j$re_rezi_a_nc)

ggsave("output/props.svg", height = 5, width = 6)