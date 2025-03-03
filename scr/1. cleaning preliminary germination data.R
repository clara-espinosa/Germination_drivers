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

# cold_stratification pretreatment geermination and mean values ####
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
  
# DARKNESS treatment special calculations to remove cold stratification germ ####
# calculate mean cold stratification germination to apply to darkness treatment
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, petri, initial_pre, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code, petri) %>%
  dplyr::summarize(cold_strat = round(mean(cold_strat),0)) %>% # mean cold stratification germination rounded without decimals
  mutate (treatment = "B_alternate_dark") %>%
  as.data.frame()-> mean_coldstrat_germ

# final germination calculation
read.csv("data/raw_data.csv", sep = ",") %>%
  filter(treatment=="B_alternate_dark")%>%
  convert_as_factor(species, code, treatment, temp)%>%
  dplyr::select(species, code, treatment, petri, initial_pre, D42, empty, fungus)%>%
  merge(mean_coldstrat_germ, by = c("species", "code", "petri", "treatment"))%>%
  group_by(species, code, treatment, petri)%>%
  summarise(viable=initial_pre-empty-fungus-cold_strat, # remove cold stratification germination
            finalgerm=D42-cold_strat)%>% # final day germination subtracting germination happened under cold stratification
  mutate (finalgerm = ifelse(finalgerm < 0, 0, finalgerm),# change negative values from calculated mean to 0 
          viable = ifelse(viable < 0, 0, viable))-> dark_germ  
  
  
# combine above datasets to have 1 working dataframe (without the cold stratification treatment) ####
read.csv("data/raw_data.csv", sep = ",") %>%
  filter(!treatment=="B_alternate_dark")%>%
  convert_as_factor(species, code, treatment, temp) %>%
  mutate (initial_treat = initial_pre-cold_strat)%>%
  mutate (viable = initial_treat -(empty + fungus)) %>%
  dplyr::select (species, code, treatment, petri, initial_treat, D7:D42,viable)%>%
  gather ("scores", "germ", D7:D42)%>%    
  group_by (species, code, treatment, petri)%>%
  dplyr::summarize (finalgerm = sum(germ, na.rm=T), 
                    viable = mean(viable, na.rm=T)) %>% 
  rbind(cold_strat)%>% # for germination under cold stratification
  rbind(dark_germ)%>%  # for darkness treatment
  mutate(treatment = fct_relevel(treatment, "A_alternate_light", "B_alternate_dark", "C_alternate_WP", "D_constant_light","E_cold_stratification" ))%>%
  arrange(species, code, treatment)-> finalgerm 

finalgerm%>%
  #filter(!treatment=="E_cold_stratification")%>%
  group_by(species, code)%>%
  summarise(finalgerm = sum(finalgerm))%>%
  #print (n=58)
  filter(!finalgerm>0)%>%
  select(species)-> nogerm_species


# 7 species (all temperate) with 0 germination across all experiment (remove from analysis)
# temperate species 
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%

finalgerm%>%
    filter(!treatment=="E_cold_stratification")%>%
    group_by(species, code)%>%
    summarise(finalgerm = sum(finalgerm))%>%
    #print (n=58)
    filter(!finalgerm>0)%>%
    select(species) 
# 4 species (all mediterranean) only germinated in cold stratification
  filter (!species == "Teesdalia conferta") # all germinated in cold stratification
  filter (!species == "Solidago virgaurea") # only germinated in cold stratification
  filter (!species == "Phalacrocarpum oppositifolium") # only germinated in cold stratification >50%
  filter (!species == "Avenella flexuosa") # only germinated in cold stratification >50%