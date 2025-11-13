library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(binom)

# dataframe with species data ####
read.csv("data/species.csv", sep=",") %>%
  convert_as_factor(species, opt_temp, family, community)-> species
str(species)
unique(species$species)

# create cold_stratification germination dataframe####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, petri, initial_pre, cold_strat) %>% ## select variables
  na.omit() %>% #remove alternate dark treatment rows (no measurement of cold stratification)
  group_by(species, code, petri) %>% # group data per species, code and Petri dish
  dplyr::summarize(finalgerm = sum(cold_strat), # new summarised variables
                   viable = sum(initial_pre)) %>% 
  mutate (finalgerm = replace_na(finalgerm , 0)) %>% 
  # salix and solidago species have some petridish with all fungus
  mutate (treatment = "E_cold_stratification") %>% # create a new variable with the treatment name
  mutate(treatment = as.factor(treatment))-> cold_strat  # new cold stratification dataframe
unique(cold_strat$species)

cold_strat%>%
  group_by(species)%>% # grouped data per species
  summarise(finalgerm=sum(finalgerm), # summarise variables
            viable=sum(viable))%>%
  mutate(germpro=finalgerm/viable)%>% # calculate the germination proportion
  filter(germpro>0.5)%>% # keep only species with above 50% germination during cold stratification
  select(species)->coldstrat_highgerm # dataframe with species >50% germ in cold strat
  
# DARKNESS treatment germination calculations removing cold stratification germ ####
# calculate mean cold stratification germination to apply to darkness treatment
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  dplyr::select (species, code, petri, initial_pre, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code, petri) %>%
  # mean cold stratification germination per species in cold strat that 
  # will be substracted to final germination in darkness
  dplyr::summarize(cold_strat = round(mean(cold_strat),0)) %>% 
  mutate (treatment = "B_alternate_dark") %>% # new variable to correctly merge dataframes
  as.data.frame()-> mean_coldstrat_germ

# final germination calculation in darkness 
read.csv("data/raw_data.csv", sep = ",") %>%
  filter(treatment=="B_alternate_dark")%>% # keep only data from darkness treatment
  convert_as_factor(species, code, treatment, temp)%>%
  dplyr::select(species, code, treatment, petri, initial_pre, D42, empty, fungus)%>%
  # merge dataframe with calculated mean germination during cold stratification per species
  merge(mean_coldstrat_germ, by = c("species", "code", "petri", "treatment"))%>% 
  group_by(species, code, treatment, petri)%>%
  # remove calculated cold stratification germination to both viable seeds and final germination recorded
  summarise(viable=initial_pre-empty-fungus-cold_strat, 
            finalgerm=D42-cold_strat)%>% # D42 day with final germination registered
  # due to using averages in cold stratification, some values might end up slightly negative (-1 or -2)
  # thus we change negative values from calculated mean to 0 germination
  mutate (finalgerm = ifelse(finalgerm < 0, 0, finalgerm),
          viable = ifelse(viable < 0, 0, viable))-> dark_germ # dataframe with final germination from the darkness treatment
  
# combine above germination datasets to have 1 working dataframe ####
read.csv("data/raw_data.csv", sep = ",") %>%
  filter(!treatment=="B_alternate_dark")%>% # remove raw dark germination data to join the proper dataset afterwards
  convert_as_factor(species, code, treatment, temp) %>%
  # remove cold stratification germination in the rest of the treatments control, water limitation and constant temperatures
  mutate (initial_treat = initial_pre-cold_strat)%>% 
  mutate (viable = initial_treat -(empty + fungus)) %>% # viable seeds calculations
  dplyr::select (species, code, treatment, petri, initial_treat, D7:D42,viable)%>%
  gather ("scores", "germ", D7:D42)%>%    # transform from wide to long format
  group_by (species, code, treatment, petri)%>%
  dplyr::summarize (finalgerm = sum(germ, na.rm=T), # sum germination per petri 
                    viable = mean(viable, na.rm=T)) %>% # keep the value of viable seeds x petri (no need to sum)
  rbind(cold_strat)%>% # dataframe for germination under cold stratification
  rbind(dark_germ)%>%  # dataframe for darkness treatment
  # arrange the treatments for easier comparison/representation
  mutate(treatment = fct_relevel(treatment, "A_alternate_light", "B_alternate_dark", "C_alternate_WP", "D_constant_light","E_cold_stratification" ))%>%
  arrange(species, code, treatment)-> finalgerm 

# "A_alternate_light" = Control treatment
# "B_alternate_dark" = Darkness treatment
# "C_alternate_WP" = Water limitation treatment
# "D_constant_light" = Constant temperatures treatment
# "E_cold_stratification" = Cold stratification 

# 7 species (all from temperate community) with 0 germination across all experiment (removed from analysis)
 finalgerm%>%
  group_by(species, code)%>%
  summarise(finalgerm = sum(finalgerm))%>%
  filter(!finalgerm>0)%>%
  select(species)-> nogerm_species # to remove in further analysis


# 10 species >50% germ in cold stratification that will also be removed from further analysis
finalgerm%>%
    filter(treatment=="E_cold_stratification")%>%
    group_by(species, code)%>%
    summarise(finalgerm = sum(finalgerm), viable=sum(viable))%>%
    mutate(germPRO=finalgerm/viable)%>%
    filter(germPRO>0.5)%>%
    select(species) -> coldstrat_highgerm # to remove in further analysis
  
# Species with more than 10% germ in control conditions
  finalgerm%>%
    filter(treatment=="A_alternate_light")%>%
    filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
    filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
    group_by(species, code)%>%
    summarise(finalgerm= sum(finalgerm), viable=sum(viable))%>%
    mutate(germPRO=finalgerm/viable)%>%
    filter(germPRO>0.1)%>%
    # change names only for appropiate phylogenetic random factor in MCMC GLMM
    mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>% 
    select(species)-> above10_control
unique(above10_control$species)

# Species with more than 50% germ in control conditions (same scrpit lines as above)
  finalgerm%>%
    filter(treatment=="A_alternate_light")%>%
    filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
    filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
    group_by(species, code)%>%
    summarise(finalgerm= sum(finalgerm), viable=sum(viable))%>%
    mutate(germPRO=finalgerm/viable)%>%
    filter(germPRO>0.5)%>%
    mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
    select(species)-> above50_control
  