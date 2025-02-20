library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(binom)

# dataframe with raw data ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) -> raw_df
str(raw_df)
unique(raw_df$species)

# dataframe with species data ####
read.csv("data/species.csv", sep=",") %>%
  dplyr::select(species, opt_temp, family, community, biogeography,habitat)  %>%
  convert_as_factor(species, opt_temp, family, community, biogeography,habitat)  -> species
str(species)
unique(species$species)
setdiff(species$species, raw_df$species)
setdiff(raw_df$species,species$species)

# D0 germination check change to cold_stratification treatment and mean values ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, petri, initial_pre, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code,petri) %>%
  dplyr::summarize(finalgerm = sum(cold_strat),
                   viable = sum(initial_pre)) %>% 
  mutate (finalgerm = replace_na(finalgerm , 0)) %>% # salix and solidago have some petridish with all fungus
  mutate (treatment = "E_cold_stratification") %>%
  mutate(treatment = as.factor(treatment))-> cold_strat 
unique(cold_strat$species)

cold_strat%>%
  group_by(species)%>%
  summarise(finalgerm=sum(finalgerm),
            viable=sum(viable))%>%
  mutate(germpro=finalgerm/viable)%>%
  filter(germpro>0.5)%>%
  select(species)->coldstrat_highgerm
  
# calculate mean cold stratification germination to apply to darkness treatment
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, petri, initial_pre, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code, petri) %>%
  dplyr::summarize(cold_strat = round(mean(cold_strat),0)) %>% # mean cold stratification germination rounded without decimals
  mutate (treatment = "B_alternate_dark") %>%
  as.data.frame()-> mean_coldstrat_germ

# combine above dataset to have 1 working dataframe (without the cold stratification treatment) ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, treatment,petri, cold_strat) %>%
  na.omit()%>%
  rbind(mean_coldstrat_germ)%>% ## ad mean germ under cold stratification per petri dish to fill darkness cold strat values
  arrange(species, treatment)%>%
  cbind(raw_df%>%select(initial_pre, D7:fungus)) %>%
  mutate (initial_treat = initial_pre-cold_strat)%>%
  mutate (viable = initial_treat -(empty + fungus)) %>%
  mutate (viable = ifelse(viable < 0, 0, viable)) %>% # change negative values from caluclated mean to 0 
  dplyr::select (species, code, treatment, petri, initial_pre, cold_strat, initial_treat, D7:D42,viable)-> germ_df 
unique(germ_df$species)

# viables x community x treatment/species
germ_df %>%
  merge(species)%>%
  group_by(community)%>%
  summarise(viable= sum(viable), 
            cold_strat=sum(cold_strat))%>%
  mutate(total_viable = viable+cold_strat)

# list of accession with less than 25% of viable seeds only salix (remove from analysis) 
germ_df %>%
  group_by (species, code, treatment)%>%
  mutate (viablePER = (viable/initial_treat)*100) %>% 
  group_by(species, code) %>%
  summarise(viable = sum(viable), viablePER = mean(viablePER))%>%
  filter(viablePER<25) # phalacrocarpum oppositifolium

# final germination percentage calculation discounting cold stratification germ ####
germ_df%>%
  convert_as_factor(species, code, treatment) %>%
  dplyr::select (species, code, treatment, petri, initial_treat, D7, D14, D21, D28, D35, D42, viable) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species, code, treatment, petri)%>%
  dplyr::summarize (finalgerm = sum(germ), 
                    viable = mean(viable)) %>% 
  rbind(cold_strat)%>%
  mutate (finalgerm = ifelse(finalgerm > viable, viable, finalgerm))%>%  # for darkness treatment
  arrange(species, code, treatment)-> finalgerm 

unique(finalgerm$species)

finalgerm%>%
  #filter(!treatment=="E_cold_stratification")%>%
  group_by(species, code)%>%
  summarise(finalgerm = sum(finalgerm))%>%
  #print (n=58)
  filter(!finalgerm>0)%>%
  select(species)-> nogerm_species
  

# species with 0 germination across all experiment (remove from analysis)
# temperate species 
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%

# Mediterranean species  
  filter (!species == "Teesdalia conferta") # all germinated in cold stratification
  filter (!species == "Solidago virgaurea") # only germinated in cold stratification
  filter (!species == "Phalacrocarpum oppositifolium") # only germinated in cold stratification
  filter (!species == "Avenella flexuosa") # only germinated in cold stratification