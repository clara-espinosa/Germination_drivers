library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(vegan); library(ade4); library(MASS)
library(ecodist);library(maptools);library(rpart);library(splines)
library(gam);library(pgirmess);library(utils);library(combinat)
library(cluster);library(fpc);library(clusterSim);library(lmtest)
library(Hmisc);library(gplots);library(NbClust);library(rpart)
library(rpart.plot);library(dismo);library(multcomp);library(gbm)
library(raster);library(tidyverse);
library(glmmTMB); library (DHARMa)

# following chapter 4 multivariate species level responses by Lars ###
# Separate from the beggining between communities
############################################# MEDITERRANEAN ###########################################################
# 4.1 DATA MATRICES ####
# header species data 
read.csv("data/species.csv") %>%
  dplyr::select(species,code, family, community, habitat, germ_drivers)%>%
  filter(community == "Mediterranean") -> spe_med

# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
spe_med$sp <- make.cepnames(spe_med$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x community matrix
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  merge(spe_med, by="species")%>% # get species names shortened
  dplyr::select(plot, sp, rel_cover)%>%
  spread(sp, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  filter(!plot=="A00")%>% # filter plots without iButtons data
  filter(!plot=="B00")%>%
  filter(!plot=="C00")%>%
  filter(!plot=="D00")%>%
  filter(!plot=="B13")%>%
  filter(!plot=="B14")%>%
  filter(!plot=="B15")%>%
  filter(!plot=="C04")%>%
  filter(!plot=="C07")%>%
  column_to_rownames(var="plot")-> sp_x_com_M # replace Na with 0

# 2-  species x traits matrix with germination proportion
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
  #filter (!species == "Solidago virgaurea")%>%  #filter species with 0 germination traits?
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  #filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  rbind(cold_strat_M)%>%
  mutate(germpro = round(finalgerm/viable, 2))%>%
  merge(spe_med, by =c("species", "code"))%>% 
  dplyr::select(sp, treatment, germpro)%>%
  spread(treatment, germpro)->germ_trait_M  
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
  #filter (!species == "Solidago virgaurea")%>%  #filter species with 0 germination traits?
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  #filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  rbind(cold_strat_M)%>%
  merge(spe_med)-> germ_to_odds_M

library(questionr)
germ_to_odds_M%>%
  filter(species=="Dianthus langeanus")-> D_test
fisher.test(D_test)

glm(cbind(finalgerm, viable - finalgerm) ~ treatment,  family = "binomial", data= D_test)-> D_odds
summary(D_odds)
odds.ratio(D_odds)
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
  merge(spe_med, by= "species") -> seed_production_M

# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  merge(spe_med, by= "species") -> seed_mass_M 

setdiff(species$species, test1$species)
# plant and floral height (only Teesdalia missing)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2),
                   floral_height = round(mean(floral_height),2))%>%
  right_join(spe_med, by= "species")-> plant_height_M 
setdiff(species$species, test1$species)
# join traits subset
seed_production_M%>%
  merge(seed_mass_M)%>%
  merge(plant_height_M)%>%
  merge(germ_odds_M)%>%
  dplyr::select(sp, seed_mass,seed_production,plant_height,floral_height, 
                odds_A_control, odds_B_dark, odds_C_WP, odds_D_constant,odds_E_cold )%>%
  replace(is.na(.), 10)%>%
  merge(germ_trait_M, by = "sp")%>%
  column_to_rownames(var="sp")->sp_x_trait_M
  
# 3- Community x Environmental data 
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
  column_to_rownames(var="plot")-> com_x_env_M
# 4.3 Trait free CCA ####
cca1M <- cca(sp_x_com_M ~ elevation+GDD+FDD, data=com_x_env_M) # +bio1+bio2+bio7+Snw tried but highly correlated
plot(cca1M, display = c("species", "bp"))
anova(cca1M, by= "terms") # elevation and GDD significantly explain the ordination of the species
RsquareAdj(cca1M) # roughly 6% of variation in the species data is explained by these variables
head(vegan::scores(cca1M, display = "species"))
# possibility to separate CCA for each env. variable and then use the first CCA axes to obtain 
# the value (score) for the particular gradient

# 4.4 Double canonical correspondance analysis (dbrda function not working!!) ####
ca1M <- dudi.coa(sp_x_com_M, scannf = F)
dCCA1M <- ade4::dbrda(ca1M, com_x_env_M, sp_x_trait_M, scannf = FALSE)
# 4.5 Functional response groups ####
dist_CCAM <- dist(vegan::scores(cca1M, display = "species"))
clust_CCAM <- hclust(dist_CCAM, method = "ward.D2")
plot(clust_CCAM, cex = 0.6)

#The function NbClust from the package of the same name estimates the optimal number of groups on the basis of maximizing 
#the dissimilarity between groups and minimizing the dissimilarity within groups.
groups_CCAM <- NbClust(diss = dist_CCAM, distance = NULL, min.nc = 2, max.nc = 6,
                       method = "ward.D2", index = "silhouette")
groups_CCAM # optimal 3 clusters

# we are asking whether there are traits that explain to which cluster the species belong, keeping in
# mind that this clustering is based on our initial ordination that arranges species along the
# considered environmental gradients. For this test we use here a simple regression analysis
# test every trait against the cluster affiliation
summary(lm(as.matrix(sp_x_trait_M) ~ groups_CCAM$Best.partition))
# No trait significantly associated to each of the groups

# 4.6 RDA and regression trees (more appropiate with short gradients and sp abundance-gradient linear relationship) ####
rda1M <- rda(sp_x_com_M ~ ., data = com_x_env_M, scale = TRUE)
plot(rda1M, display = c("bp", "sp"))
RsquareAdj(rda1M)
#how much of the variation the 3 constraining factors explain.
cumsum(rda1M$CCA$eig) / sum(rda1M$CCA$eig)

# Analyzing the relationship between the species responses (expressed as the species scores in
#the RDA ordination space) and their traits
lm_rda1M <- lm(vegan::scores(rda1M, choices = 1, display = "species") ~ ., data = sp_x_trait_M)
summary(lm_rda1M)
# significance of odds_D_constant
# marginal significance seed mass, floral height, D_constant_proportion
# the above traits of our species do no explain any significant proportion of the sp responses
### STOPED FOLLOWING CHAPTER 4 HERE, NEXT WOULD BE REGRESSION TREES

################################# TEMPERATE ########################################################################
# 4.1 DATA MATRICES ####
# header species data 
read.csv("data/species.csv") %>%
  dplyr::select(species, code, family, community, habitat, germ_drivers)%>%
  filter(community == "Temperate")%>%
  filter (!species == "Minuartia CF")-> spe_tem
# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
spe_tem$sp <- make.cepnames(spe_tem$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x community matrix
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species))) %>%
  merge(spe_tem, by="species")%>% # get species names shortened
  dplyr::select(plot, sp, rel_cover)%>%
  spread(sp, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>%
  column_to_rownames(var="plot")-> sp_x_com_T # replace Na with 0
setdiff(spe_tem$species, test1$species)
# 2-  species x traits matrix
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
  #filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  #filter (!species == "Gentiana verna")%>%
  #filter (!species == "Gentianella campestris")%>%
  #filter (!species == "Kobresia myosuroides")%>%
  #filter (!species == "Salix breviserrata")%>%
  #filter (!species == "Sedum album")%>%
  #filter (!species == "Sedum atratum")%>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  rbind(cold_strat_T)%>%
  mutate(germpro = round(finalgerm/viable, 2))%>%
  merge(spe_tem, by =c("species", "code"))%>% 
  dplyr::select(sp, treatment, germpro)%>%
  spread(treatment, germpro)->germ_trait_T
  column_to_rownames(var="sp")%>%   #sp_x_trait_T   # sp_x_trait_M
    
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
  #filter (!species == "Solidago virgaurea")%>%  #filter species with 0 germination traits?
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
  right_join(spe_tem, by= "species") -> seed_production_T
setdiff(spe_tem$species, test1$species)
# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  merge(spe_tem, by= "species") -> seed_mass_T 
  
setdiff(spe_tem$species, test1$species)
# plant and floral height (only Gentianella campestris missing)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2),
                     floral_height = round(mean(floral_height),2))%>%
  right_join(spe_tem, by= "species")-> plant_height_T
setdiff(spe_tem$species, test1$species)

# join traits subset
seed_production_T%>%
  merge(seed_mass_T)%>%
  merge(plant_height_T)%>%
  merge(germ_odds_T)%>%
  dplyr::select(sp, seed_mass,seed_production,plant_height,floral_height, 
                odds_A_control, odds_B_dark, odds_C_WP, odds_D_constant,odds_E_cold )%>%
  replace(is.na(.), 10)%>% # to fill Na until complete all traits
  merge(germ_trait_T, by = "sp")%>%
  column_to_rownames(var="sp")->sp_x_trait_T

# 3- Community x Environmental data
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
  column_to_rownames(var="plot")-> com_x_env_T

# 4.3 Trait free CCA ####
cca1T <- cca(sp_x_com_T ~ elevation+GDD+FDD+bio1+bio2+bio7+Snw, data=com_x_env_T) #  tried but highly correlated
plot(cca1T, display = c("species", "bp"))
anova(cca1T, by= "terms") # elevation, GDD, FDD, bio7 and snow significantly explain the ordination of the species
RsquareAdj(cca1T) # roughly 15% of variation in the species data is explained by these variables
head(vegan::scores(cca1T, display = "species"))
# separate CCA for each env. variable and then use the first CCA axes to obtain the value for the particular gradient
# still needs to account for the effects of the other variables as conditional variables (covariates)

# 4.4 Double canonical correspondance analysis (dbrda function not working!!) ####
ca1T<- dudi.coa(sp_x_com_T, scannf = F)
dCCA1T <- ade4::dbrda(ca1T, com_x_env_T, sp_x_trait_T, scannf = FALSE)

# 4.5 Functional response groups ####
dist_CCAT <- dist(vegan::scores(cca1T, display = "species"))
clust_CCAT <- hclust(dist_CCAT, method = "ward.D2")
plot(clust_CCAT, cex = 0.6)

#The function NbClust from the package of the same name estimates the optimal number of groups on the basis of maximizing 
#the dissimilarity between groups and minimizing the dissimilarity within groups.
groups_CCAT <- NbClust(diss = dist_CCAT, distance = NULL, min.nc = 2, max.nc = 6,
                       method = "ward.D2", index = "silhouette")
groups_CCAT# optimal 2 clusters

# we are asking whether there are traits that explain to which cluster the species belong, keeping in
# mind that this clustering is based on our initial ordination that arranges species along the
# considered environmental gradients. For this test we use here a simple regression analysis
# test every trait against the cluster affiliation
summary(lm(as.matrix(sp_x_trait_T) ~ groups_CCAT$Best.partition))
# No seed trait significantly associated to each of the groups

# 4.6 RDA and regression trees (more appropiate with short gradients and sp abundance-gradient linear relationship) ####
rda1T <- rda(sp_x_com_T ~ ., data = com_x_env_T, scale = TRUE)
plot(rda1T, display = c("bp", "sp"))
RsquareAdj(rda1T)
#how much of the variation the 3 constraining factors explain.
cumsum(rda1T$CCA$eig) / sum(rda1T$CCA$eig)

# Analyzing the relationship between the species responses (expressed as the species scores in
#the RDA ordination space) and their traits
lm_rda1T <- lm(vegan::scores(rda1T, choices = 1, display = "species") ~ ., data = sp_x_trait_T)
summary(lm_rda1T)
# significance of seed mass, plant height, odds_A_control