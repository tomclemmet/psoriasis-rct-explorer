library(ggplot2)

pasi |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |> 
  select(ref_id, drug, pasi50, pasi75, pasi90, pasi100) |> 
  group_by(ref_id) |> 
  summarise(
    p50 = any(!is.na(pasi50) & is.na(pasi75) & is.na(pasi90) & is.na(pasi100)),
    p75 = any(is.na(pasi50) & !is.na(pasi75) & is.na(pasi90) & is.na(pasi100)),
    p90 = any(is.na(pasi50) & is.na(pasi75) & !is.na(pasi90) & is.na(pasi100)),
    p100 = any(is.na(pasi50) & is.na(pasi75) & is.na(pasi90) & !is.na(pasi100)),
    p50_75 = any(!is.na(pasi50) & !is.na(pasi75) & is.na(pasi90) & is.na(pasi100)),
    p50_90 = any(!is.na(pasi50) & is.na(pasi75) & !is.na(pasi90) & is.na(pasi100)),
    p50_100 = any(!is.na(pasi50) & !is.na(pasi75) & is.na(pasi90) & !is.na(pasi100)),
    p75_90 = any(is.na(pasi50) & !is.na(pasi75) & !is.na(pasi90) & is.na(pasi100)),
    p75_100 = any(is.na(pasi50) & !is.na(pasi75) & is.na(pasi90) & !is.na(pasi100)),
    p90_100 = any(is.na(pasi50) & is.na(pasi75) & !is.na(pasi90) & !is.na(pasi100)),
    p50_75_90 = any(!is.na(pasi50) & !is.na(pasi75) & !is.na(pasi90) & is.na(pasi100)),
    p50_90_100 = any(!is.na(pasi50) & is.na(pasi75) & !is.na(pasi90) & !is.na(pasi100)),
    p50_75_100 = any(!is.na(pasi50) & !is.na(pasi75) & is.na(pasi90) & !is.na(pasi100)),
    p75_90_100 = any(is.na(pasi50) & !is.na(pasi75) & !is.na(pasi90) & !is.na(pasi100)),
    p50_75_90_100 = any(!is.na(pasi50) & !is.na(pasi75) & !is.na(pasi90) & !is.na(pasi100))
  ) |> 
  ungroup() |> 
  summarise(across(p50:p50_75_90_100, sum))


pasi |> 
  filter(!if_all(pasi50:pasi100, \(x) is.na(x))) |> 
  select(ref_id, drug, n, pasi50, pasi75, pasi90, pasi100) |> 
  pivot_longer(pasi50:pasi100, names_to = "outcome", values_to = "r") |> 
  filter(!is.na(r)) |> 
  summarise(.by = c(ref_id, outcome, drug), n = sum(n), r = sum(r)) |> 
  mutate(p = r / n) |> 
  ggplot(aes(x = p, y = drug)) +
  geom_point(aes(size = n, colour = outcome), alpha = 0.5) +
  facet_wrap(~ outcome)
