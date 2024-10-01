library(ggplot2);library(tidyverse);library(rstatix);library(lubridate);library (vegan)
# Supporting information

#2-  species x traits  MEDITERRANEAN matrix ########
#with germination proportion (NOT TO BE USED, REMOVE?)
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, treatment, petri, initial, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species, code, treatment, petri)%>%
  dplyr::summarize (finalgerm = sum(germ)) %>% # 
  dplyr::select(species, code, treatment, petri, finalgerm)%>% 
  merge(viables_petri, by= c("species", "code", "treatment", "petri")) %>%
  rbind(cold_strat)%>%
  group_by (species, code, treatment) %>%
  dplyr::summarize(finalgerm = sum (finalgerm)) %>% 
  filter (!species == "Solidago virgaurea")%>%  #filter species with 0 germination traits?
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  #filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  rbind(cold_strat_M)%>%
  mutate(germper = round((finalgerm/viable)*100, 2))%>%
  mutate(sqrt_germper = sqrt(germper))%>%
  merge(spe_med, by =c("species", "code"))%>% 
  dplyr::select(species, treatment, sqrt_germper)%>%
  spread(treatment, sqrt_germper)%>%
  rename(A_control = A_alternate_light)%>%
  rename(B_dark = B_alternate_dark)%>%
  rename(C_WP = C_alternate_WP)%>%
  rename(D_constant = D_constant_light)%>%
  rename(E_cold = E_cold_stratification)-> germ_trait_M  
setdiff(spe_med$species,test1$species)
unique(germ_trait_M  $species)
#2-  species x traits  TEMPERATE matrix ########
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, treatment, petri, initial, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species, code, treatment, petri)%>%
  dplyr::summarize (finalgerm = sum(germ)) %>% # 
  dplyr::select(species, code, treatment, petri, finalgerm)%>% 
  merge(viables_petri, by= c("species", "code", "treatment", "petri")) %>%
  rbind(cold_strat)%>%
  group_by (species, code, treatment) %>%
  dplyr::summarize(finalgerm = sum (finalgerm)) %>%  
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  rbind(cold_strat_T)%>%
  mutate(germper = round((finalgerm/viable)*100, 2))%>%
  mutate(sqrt_germper = sqrt(germper))%>%
  merge(spe_tem, by =c("species", "code"))%>% 
  dplyr::select(species, treatment, sqrt_germper)%>%
  spread(treatment, sqrt_germper)%>%
  rename(A_control = A_alternate_light)%>%
  rename(B_dark = B_alternate_dark)%>%
  rename(C_WP = C_alternate_WP)%>%
  rename(D_constant = D_constant_light)%>%
  rename(E_cold = E_cold_stratification)->germ_trait_T