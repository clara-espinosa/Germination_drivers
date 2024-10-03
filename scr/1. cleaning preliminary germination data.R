library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(binom)

# dataframe with raw data ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) -> raw_df
str(raw_df)
unique(raw_df$species)

# dataframe with species data ####
read.csv("data/species.csv", sep=",") %>%
  dplyr::dplyr::select(species, code, temperature_regime, family, community, habitat)  %>%
  convert_as_factor(species, code, temperature_regime, family, community, habitat)  -> species
str(species)
unique(species$species)

# species with 0 germination across all experiment (remove from analysis)
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Teesdalia conferta") # all germinated in cold stratification
  filter (!species == "Solidago virgaurea")
  
# viables calculation x petri dish ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  dplyr::select (species, code, treatment, petri, viable) -> viables_petri
unique(viables_petri$species)
# viables x community x treatment
read.csv("data/raw_data.csv", sep = ",") %>%
    convert_as_factor(species, code, treatment, temp) %>%
    dplyr::select (species, code, treatment, petri, initial, empty, fungus) %>%
    group_by (species, code, treatment, petri) %>%
    mutate (viable = initial -(empty + fungus)) %>%
    dplyr::select (species, code, treatment, petri, viable)%>%
    merge(species)%>%
    group_by (community, treatment)%>%
    summarise(viable= sum(viable))-> viables_community

# viables x species x treatment
read.csv("data/raw_data.csv", sep = ",") %>%
    convert_as_factor(species, code, treatment, temp) %>%
    dplyr::select (species, code, treatment, petri, initial, empty, fungus) %>%
    group_by (species, code, treatment, petri) %>%
    mutate (viable = initial -(empty + fungus)) %>%
    dplyr::select (species, code, treatment, petri, viable)%>%
    group_by (species,code, treatment)%>%
    summarise(viable= sum(viable))-> viables_sp
unique(viables_sp$species)

# list of accession with less than 25% of viable seeds (remove from analysis) 
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  mutate (viablePER = (viable/initial)*100) %>% 
  group_by(species, code, treatment) %>%
  summarise(viable = sum(viable), viablePER = mean(viablePER))%>%
  filter(viablePER<25) # salix breviserrata


# D0 germination check change to cold_stratification treatment and mean values ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  merge(viables_petri) %>%
  dplyr::select (species, code, petri, viable, initial, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code,petri) %>%
  dplyr::summarize(finalgerm = sum(cold_strat),
            viable = sum(viable)) %>% 
  mutate (finalgerm = replace_na(finalgerm , 0)) %>% # salix and solidago have some petridish with all fungus
  mutate (treatment = "E_cold_stratification") %>%
  mutate(treatment = as.factor(treatment))-> cold_strat 
unique(cold_strat$species)

#needed?
cold_strat %>%
  group_by(species, code, treatment)%>%
  dplyr::summarize(finalgerm = sum(finalgerm),
                   viable = sum(viable))%>%
  merge(species)%>%
  filter(community=="Mediterranean")%>%
  dplyr::select(code, species, treatment, viable, finalgerm)->cold_strat_M
#needed?
cold_strat %>%
  group_by(species, code, treatment)%>%
  dplyr::summarize(finalgerm = sum(finalgerm),
                   viable = sum(viable))%>%
  merge(species)%>%
  filter(community=="Temperate")%>%
  dplyr::select(code, species, treatment, viable, finalgerm)->cold_strat_T


# final germination percentage calculation discounting cold stratification germ ####
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
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))-> finalgerm 
