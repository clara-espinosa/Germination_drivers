library(ggplot2);library(tidyverse);library(rstatix);library(lubridate)
# Data matrices handling
############################################# MEDITERRANEAN ###########################################################
# DATA MATRICES ####
# header species data ####
read.csv("data/species.csv") %>%
  dplyr::select(species,code, family, community, habitat, germ_drivers)%>%
  filter(community == "Mediterranean")%>%
  filter(!species=="Solidago virgaurea")-> spe_med

# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
spe_med$sp <- make.cepnames(spe_med$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x plot matrix #####
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  filter(!species=="Solidago virgaurea")%>%
  merge(spe_med, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  spread(species, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  filter(!plot=="A00")%>% # filter plots without iButtons data or without any vegetation
  filter(!plot=="B00")%>%
  filter(!plot=="C00")%>%
  filter(!plot=="D00")%>%
  filter(!plot=="B13")%>%
  filter(!plot=="B14")%>%
  filter(!plot=="B15")%>%
  filter(!plot=="C04")%>%
  filter(!plot=="C07")%>%
  column_to_rownames(var="plot")-> plot_x_sp_M # replace Na with 0
sp_x_plot_M <- t(plot_x_sp_M)

# presence_absence sp x plot matrix
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  filter(!species=="Solidago virgaurea")%>%
  merge(spe_med, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  filter(!plot=="A00")%>% # filter plots without iButtons data or without any vegetation
  filter(!plot=="B00")%>%
  filter(!plot=="C00")%>%
  filter(!plot=="D00")%>%
  filter(!plot=="B13")%>%
  filter(!plot=="B14")%>%
  filter(!plot=="B15")%>%
  filter(!plot=="C04")%>%
  filter(!plot=="C07")%>%
  column_to_rownames(var="plot")-> plot_x_sp_M_PA 

# 2-  species x traits matrix ########
#with germination proportion
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
unique(viables_petri$species)

#### species x traits matrix with germination  as odd ratios 
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
  merge(spe_med)-> germ_to_odds_M

### Get GLM coefficients since the coefficients(estimate) are returned in log odds, 
#https://www.r-bloggers.com/2016/11/introducing-r-package-oddsratio/

glms <- function(x) {
  glm(cbind(finalgerm, viable - finalgerm) ~ treatment,  family = "binomial", data= x) -> m1
  broom::tidy(m1)
}

germ_to_odds_M %>%
  group_by(code, species) %>%
  do(glms(.)) %>%
  dplyr:: select(code, species, term, estimate) %>%
  rename(log_odds_ratio = estimate ) %>%
  mutate(term = as.factor(term))%>%
  mutate(term = fct_recode(term, odds_A_control = "(Intercept)", odds_B_dark = "treatmentB_alternate_dark",
                           odds_C_WP = "treatmentC_alternate_WP", odds_D_constant = "treatmentD_constant_light",
                           odds_E_cold = "treatmentE_cold_stratification"))%>%
  dplyr::group_by(species)%>%
  tidyr::spread(term, log_odds_ratio)%>%
  merge(spe_med, by= c("species", "code"))-> germ_odds_M

# seed production trait (all sp covered)
read.csv("data/seed_production.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(seed_production= round(mean(N.seeds),2))%>%
  mutate(seed_production = log(seed_production))%>%
  merge(spe_med, by= "species") -> seed_production_M
hist(seed_production_M$seed_production)

# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  mutate(seed_mass= log (seed_mass))%>%
  merge(spe_med, by= "species") -> seed_mass_M 
hist(seed_mass_M$seed_mass)
setdiff(species$species, test1$species)
# plant and floral height (only Teesdalia missing)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2),
                   floral_height = round(mean(floral_height),2))%>%
  mutate(plant_height=log(plant_height),
         floral_height = log (floral_height))%>%
  right_join(spe_med, by= "species")-> plant_height_M 
hist(plant_height_M$plant_height)
hist(plant_height_M$floral_height)
setdiff(species$species, test1$species)
# join traits subset
germ_trait_M%>% # sqrt of germination percetanges
  merge(germ_odds_M)%>% # log odds ratios
  merge(seed_mass_M)%>% # log transformed
  merge(plant_height_M)%>% # log transformed
  merge(seed_production_M)%>% #log transformed
  dplyr::select(species, seed_mass,seed_production,plant_height,floral_height, 
                odds_B_dark, odds_C_WP, odds_D_constant,odds_E_cold, # odds_A_control not necesary
                A_control,B_dark,C_WP,D_constant,E_cold)%>%
  column_to_rownames(var="species")->sp_x_trait_M

# 3-  plot x Environmental data ########
read.csv("data/spatial-survey-temperatures-Med.csv") %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M")) %>%
  group_by(plot, Day = lubridate::floor_date(Time, "day")) %>%
  dplyr::summarize(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 5, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  group_by(plot, Month = lubridate::floor_date(Day, "month")) %>%
  dplyr::summarize(T = mean(T), X = mean(X), N = mean(N), # Daily mean, max, min
                   Snow = sum(Snow), # Snow days per month
                   FreezeThaw = sum(FreezeThaw), # Freeze-thaw days per month
                   FDD = sum(FDD), # FDD per month
                   GDD = sum(GDD)) %>% # GDD per month
  group_by(plot) %>%
  dplyr::summarize(bio1 = round(mean(T),2), # Annual Mean Temperature
                   bio2 = round(mean(X - N),2), # Mean Diurnal Range (Mean of monthly (max temp - min temp))
                   bio7 = round(max(X) - min(N), 2), # Temperature Annual Range (BIO5-BIO6)
                   Snw = sum(Snow),
                   FDD = round(abs(sum(FDD)),2), # FDD per year
                   GDD = round(sum(GDD)),2) %>%  # GDD per year
  merge(read.csv("data/spatial-survey-header-Med.csv")) %>% 
  #dplyr::select(plot, site, aspect, elevation, latitude, longitude, bio1:GDD)%>% 
  dplyr::select(plot, elevation,bio1:GDD)%>% 
  filter(!plot=="D06")%>% # filter plots without vegetation
  filter(!plot=="D07")%>%
  column_to_rownames(var="plot")-> plot_x_env_M

################################# TEMPERATE ########################################################################
# DATA MATRICES ####
# header species data  #######
read.csv("data/species.csv") %>%
  dplyr::select(species, code, family, community, habitat, germ_drivers)%>%
  filter(community == "Temperate")%>%
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Minuartia CF")-> spe_tem
# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
spe_tem$sp <- make.cepnames(spe_tem$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x plot matrix #########
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species))) %>%
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Minuartia CF")%>%
  merge(spe_tem, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  spread(species, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>%
  column_to_rownames(var="plot")-> plot_x_sp_T # replace Na with 0
sp_x_plot_T <- t(plot_x_sp_T)
setdiff(spe_tem$species, test1$species)

# presence_absence sp x plot matrix
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Minuartia CF")%>%
  merge(spe_tem, by="species")%>%  # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  column_to_rownames(var="plot")-> plot_x_sp_T_PA 
# 2-  species x traits matrix ###########
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
#column_to_rownames(var="species")%>%   

#### species x traits matrix with germination  as odd ratios 
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
  #filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  rbind(cold_strat_T)%>%
  merge(spe_tem)-> germ_to_odds_T

### Get GLM coefficients since the coefficients(estimate) are returned in log odds,
#https://www.r-bloggers.com/2016/11/introducing-r-package-oddsratio/
glms <- function(x) {
  glm(cbind(finalgerm, viable - finalgerm) ~ treatment,  family = "binomial", data= x) -> m1
  broom::tidy(m1)
}

germ_to_odds_T %>%
  group_by(code, species) %>%
  do(glms(.)) %>%
  dplyr:: select(code, species, term, estimate) %>%
  rename(log_odds_ratio = estimate ) %>%
  mutate(term = as.factor(term))%>%
  mutate(term = fct_recode(term, odds_A_control = "(Intercept)", odds_B_dark = "treatmentB_alternate_dark",
                           odds_C_WP = "treatmentC_alternate_WP", odds_D_constant = "treatmentD_constant_light",
                           odds_E_cold = "treatmentE_cold_stratification"))%>%
  dplyr::group_by(species)%>%
  tidyr::spread(term, log_odds_ratio)%>%
  merge(spe_tem, by= c("species", "code"))-> germ_odds_T


# seed production trait (all sp covered except galium and salix)
read.csv("data/seed_production.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(seed_production= round(mean(N.seeds),2))%>%
  mutate(seed_production=log(seed_production))%>%
  right_join(spe_tem, by= "species") -> seed_production_T
hist(seed_production_T$seed_production)
setdiff(spe_tem$species, test1$species)
# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  mutate(seed_mass=log(seed_mass))%>%
  merge(spe_tem, by= "species") -> seed_mass_T 
hist(seed_mass_T $seed_mass) 
setdiff(spe_tem$species, test1$species)
# plant and floral height (only Gentianella campestris missing)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2),
                   floral_height = round(mean(floral_height),2))%>%
  mutate(plant_height=log(plant_height),
         floral_height=log(floral_height))%>%
  right_join(spe_tem, by= "species")-> plant_height_T
setdiff(spe_tem$species, test1$species)
hist(plant_height_T$plant_height) 
hist(plant_height_T$floral_height) 
# join traits subset
germ_trait_T%>% # sqrt germination percentages
  merge(germ_odds_T)%>% # log odd ratios
  merge(seed_mass_T)%>% #log seed mass
  merge(plant_height_T)%>% # log plant/floral height
  merge(seed_production_T)%>% # log seed_production
  dplyr::select(species, seed_mass,seed_production,plant_height,floral_height, 
                odds_B_dark, odds_C_WP, odds_D_constant,odds_E_cold, # odds_A_control not necesary
                A_control,B_dark,C_WP,D_constant,E_cold)%>%
  replace(is.na(.), 10)%>% # IMPORTANT 2 Na's seed production galium and floral height festuca rubra invested 10 value
  column_to_rownames(var="species")->sp_x_trait_T

# 3-  plot x Environmental data ##########
read.csv("data/spatial-survey-temperatures-Tem.csv") %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M")) %>%
  group_by(plot, Day = lubridate::floor_date(Time, "day")) %>%
  dplyr::summarize(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 5, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  group_by(plot, Month = lubridate::floor_date(Day, "month")) %>%
  dplyr::summarize(T = mean(T), X = mean(X), N = mean(N), # Daily mean, max, min
                   Snow = sum(Snow), # Snow days per month
                   FreezeThaw = sum(FreezeThaw), # Freeze-thaw days per month
                   FDD = sum(FDD), # FDD per month
                   GDD = sum(GDD)) %>% # GDD per month
  group_by(plot) %>%
  dplyr::summarize(bio1 = round(mean(T),2), # Annual Mean Temperature
                   bio2 = round(mean(X - N),2), # Mean Diurnal Range (Mean of monthly (max temp - min temp))
                   bio7 = round(max(X) - min(N), 2), # Temperature Annual Range (BIO5-BIO6)
                   Snw = sum(Snow),
                   FDD = round(abs(sum(FDD)),2), # FDD per year
                   GDD = round(sum(GDD)),2) %>%  # GDD per year
  merge(read.csv("data/spatial-survey-header-Tem.csv")) %>% 
  dplyr::select(plot, elevation,bio1:GDD)%>%
  #dplyr::select(plot, elevation, FDD, GDD)%>%
  filter(plot%in%tem_plot$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_env_T