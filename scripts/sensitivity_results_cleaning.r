library(tidyverse)

corn.sens.df <- readRDS("results-raw/cornsensresults.rds")
soy.sens.df <- readRDS("results-raw/soysensresults.rds")

# find peak corn
corn.twso <- corn.sens.df %>%
  pivot_longer(cols=DVS.WL:SM.PP, names_sep = "\\.",
               names_to = c("OUTSTAT","MODEL")) %>%
  filter(OUTSTAT=="TWSO") %>% rename(TWSO = "value")

corn.harvest <- corn.twso %>%
  group_by(soil, deltaT, deltaRAIN, MODEL) %>%
  mutate(harvest = !is.na(TWSO) & TWSO>0 & is.na(lead(TWSO,order_by = day))) %>%
  ungroup() %>% filter(harvest==TRUE, day > "2022-05-17") %>%
  mutate(bshl_acr = TWSO*0.039368/2.47105)


# find peak soy
soy.twso <- soy.sens.df %>%
  pivot_longer(cols=DVS.WL:SM.PP, names_sep = "\\.",
               names_to = c("OUTSTAT","MODEL")) %>%
  filter(OUTSTAT=="TWSO") %>% rename(TWSO = "value")

soy.harvest <- soy.twso %>%
  group_by(soil, deltaT, deltaRAIN, MODEL) %>%
  mutate(harvest = !is.na(TWSO) & TWSO>0 & is.na(lead(TWSO,order_by = day))) %>%
  ungroup() %>% filter(harvest==TRUE, day > "2022-05-17") %>%
  mutate(bshl_acr = TWSO*0.0367437/2.47105)

# write out results
saveRDS(corn.harvest, "results/sensitivity_corn.rds")
saveRDS(soy.harvest, "results/sensitivity_soy.rds")
