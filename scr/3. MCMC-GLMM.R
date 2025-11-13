library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(viridis);library(MCMCglmm); library(binom);library (ggpubr)
library(glmmTMB); library(DHARMa);library(scales)

### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL MCMC GLMM ####
### Read phylogeny tree
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set priors for germination models (as many prior as random factors)
priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))
  
### 1. Cold stratification differences ####
cold_strat%>% # dataframe from script 1. Handling germination data
  # merge information about species and community
  merge(read.csv("data/species.csv", sep = ";"), by = c("species", "code"))%>% 
  filter(!(species%in%nogerm_species$species))%>% # remove species with 0 germ across all experiment
  select(community, species,  treatment, petri, finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>%  # necessary to include phylogeny as random factor
  convert_as_factor(species, community, petri, collection_year)-> mcmc_cold

str(mcmc_cold) 
unique(mcmc_cold$species)# 50 species

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ community, 
                   random = ~ animal + ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_cold,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_cold2

x11()
plot(mc_cold) # diagnostics plots look good
summary(mc_cold) # model results

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_cold$VCV[,"animal"]/(mc_cold$VCV[,"animal"] + mc_cold$VCV[,"units"]) 

mean(mc_cold$VCV[,"animal"]/(mc_cold$VCV[,"animal"] + mc_cold$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_cold$VCV[,"animal"]/(mc_cold$VCV[,"animal"] + mc_cold$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_cold$VCV[,"animal"]/(mc_cold$VCV[,"animal"] + mc_cold$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(mc_cold)$Gcovariances[1, 1] %>% round(2) 
summary(mc_cold)$Gcovariances[1, 2] %>% round(2) 
summary(mc_cold)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_cold)$Gcovariances[2, 1] %>% round(2)
summary(mc_cold)$Gcovariances[2, 2] %>% round(2) 
summary(mc_cold)$Gcovariances[2, 3] %>% round(2) 

### 2. Water limitation differences ####
#2.1 Control vs water limitation germination ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep = ","), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  filter(species%in%above10_control$species)%>%
  #filter(species%in%above50_control$species)%>%
  convert_as_factor(species, community, petri, community, collection_year)%>% 
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  droplevels()-> mcmc_water
str(mcmc_water) 
unique(mcmc_water$species)# 

#descriptive stats
mcmc_water%>%
  group_by(treatment)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, #
                   random = ~ animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_water,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_water_10

x11()
plot(mc_water_10) # diagnostic plots look good
summary(mc_water_10)  # model results

plot(mc_water_50) # diagnostic plots look good 
summary(mc_water_50) # model results

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_water_10$VCV[,"animal"]/(mc_water_10$VCV[,"animal"] + mc_water_10$VCV[,"units"]) 
mean(mc_water_10$VCV[,"animal"]/(mc_water_10$VCV[,"animal"] + mc_water_10$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_water_10$VCV[,"animal"]/(mc_water_10$VCV[,"animal"] + mc_water_10$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_water_10$VCV[,"animal"]/(mc_water_10$VCV[,"animal"] + mc_water_10$VCV[,"units"]))[, 2] %>% round(2)


# Random effects animal
summary(mc_water_10)$Gcovariances[1, 1] %>% round(2) 
summary(mc_water_10)$Gcovariances[1, 2] %>% round(2) 
summary(mc_water_10)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_water_10)$Gcovariances[2, 1] %>% round(2)
summary(mc_water_10)$Gcovariances[2, 2] %>% round(2) 
summary(mc_water_10)$Gcovariances[2, 3] %>% round(2) 

#2.2 Germination reduction (%) due to water limitation compared to control ####
# Response to water stress calculated as the percentual reduction compared to control.
finalgerm %>%
  merge(read.csv("data/species.csv", sep=","), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>% # keep only the interested treatments
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  filter(species%in%above10_control$species)%>%
  #filter(species%in%above50_control$species)%>%
  mutate(germpro = finalgerm/viable)%>% # generate Nan because of 0 viable
  mutate_all(~replace(., is.nan(.), 0))%>% # change NAN to 0 for analysis
  select(community, species, habitat, community, treatment, petri,germpro, collection_year)%>%
  spread(treatment, germpro)%>%
  # calculation of germination reduction
  mutate(germ_reduction=(1-(C_alternate_WP/A_alternate_light))*100)%>%
  mutate_all(~replace(., is.nan(.), 0))%>% # change NAN to 0 for analysis
  convert_as_factor(species, community, petri, community, collection_year)%>% 
  #problem with negative values in those cases (2) were germination in water limitation was slightly higher than in control
  mutate(germ_reduction=ifelse(germ_reduction<0, 0, germ_reduction))%>% 
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% # necesary to include phylogeny as random factor
  droplevels()%>% 
  as.data.frame()-> mcmc_waterreduction2

str(mcmc_waterreduction2) # with exponential family
unique(mcmc_waterreduction2$species)
hist(mcmc_waterreduction2$germ_reduction)

#descriptive stats
mcmc_waterreduction2%>%
  group_by(species, community)%>%
  summarise(germ_reduction= mean(germ_reduction))%>%
  group_by(species)%>%
  get_summary_stats (germ_reduction)%>%
  arrange(mean)

# Priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
# TEST
MCMCglmm::MCMCglmm(germ_reduction ~ community, #
                   random = ~ animal + ID , #
                   family = "exponential", pedigree = nnls_orig, prior = priors, data = mcmc_waterreduction2,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_water_reduction2_50
x11()
plot(mc_water_reduction2_10) # diagnostics plots look good
plot(mc_water_reduction2_50) # diagnostics plots look good

summary(mc_water_reduction2_10) # model results with >10% germ in control
summary(mc_water_reduction2_50) # model results with >50% germ in control

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
lambda <- mc_water_reduction2_10$VCV[,"animal"]/(mc_water_reduction2_10$VCV[,"animal"] + mc_water_reduction2_10$VCV[,"units"]) 
mean(mc_water_reduction2_10$VCV[,"animal"]/(mc_water_reduction2_10$VCV[,"animal"] + mc_water_reduction2_10$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_water_reduction2_10$VCV[,"animal"]/(mc_water_reduction2_10$VCV[,"animal"] + mc_water_reduction2_10$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_water_reduction2_10$VCV[,"animal"]/(mc_water_reduction2_10$VCV[,"animal"] + mc_water_reduction2_10$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(mc_water_reduction2_10)$Gcovariances[1, 1] %>% round(2) 
summary(mc_water_reduction2_10)$Gcovariances[1, 2] %>% round(2) 
summary(mc_water_reduction2_10)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_water_reduction2_10)$Gcovariances[2, 1] %>% round(2)
summary(mc_water_reduction2_10)$Gcovariances[2, 2] %>% round(2) 
summary(mc_water_reduction2_10)$Gcovariances[2, 3] %>% round(2) 

### 3. Darkness and Constant temperatures cues ####
# 3.1 average responses to germination cues compared to control ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep = ","), by = c("species", "code"))%>%
  filter(!treatment== "E_cold_stratification")%>% # not necessary in this analysis
  filter(!treatment== "C_alternate_WP")%>% # not necessary in this analysis
  select(community, species, habitat, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  filter(species%in%above10_control$species)%>% # 30 species
  #filter(species%in%above50_control$species)%>% # 16 species
  convert_as_factor(species, community, petri,opt_temp, community, collection_year)%>% 
  droplevels()-> mcmc

str(mcmc) 
unique(mcmc$species)

#descriptive stats
mcmc%>%
  group_by(treatment)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

mcmc%>%
  group_by(treatment, community)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

# USE MULTINOMIAL PRIORS ####

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, #
                   random = ~ animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_treat_50

x11()
plot(mc_treat_10) # plot diagnostic look good
summary(mc_treat_10) # model results >10% germ in control

plot(mc_treat_50) # plot diagnostic look good
summary(mc_treat_50) # model results >50% germ in control

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_treat_10$VCV[,"animal"]/(mc_treat_10$VCV[,"animal"] + mc_treat_10$VCV[,"units"]) 

mean(mc_treat_10$VCV[,"animal"]/(mc_treat_10$VCV[,"animal"] + mc_treat_10$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_treat_10$VCV[,"animal"]/(mc_treat_10$VCV[,"animal"] + mc_treat_10$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_treat_10$VCV[,"animal"]/(mc_treat_10$VCV[,"animal"] + mc_treat_10$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(mc_treat_10)$Gcovariances[1, 1] %>% round(2) 
summary(mc_treat_10)$Gcovariances[1, 2] %>% round(2) 
summary(mc_treat_10)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_treat_10)$Gcovariances[2, 1] %>% round(2)
summary(mc_treat_10)$Gcovariances[2, 2] %>% round(2) 
summary(mc_treat_10)$Gcovariances[2, 3] %>% round(2) 

# 3.2 Comparison of community within each germination treatment ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep= ","), by = c("species", "code"))%>%
  # create different subsets to analyse specific responses to each treatment
  filter(treatment== "A_alternate_light" )%>% # "B_alternate_dark"  # "D_constant_light"  "A_alternate_light"
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  #filter(species%in%above10_control$species)%>%
  #filter(species%in%above50_control$species)%>%
  convert_as_factor(species, community, petri,opt_temp, community, collection_year)%>% 
  droplevels()-> mcmc_control

str(mcmc_control) 
unique(mcmc_control$species)
unique(mcmc_constant$species)
unique(mcmc_dark$species)

#descriptive
mcmc_constant%>% # mcmc_control  or mcmc_dark
  group_by(community)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

# USE MULTINOMIAL PRIORS!!!!
### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ community, #
                   random = ~ animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_control,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_control_50

x11()
plot(mc_control_10) # plot diagnostics look good
plot(mc_dark_10) # plot diagnostics look good
plot(mc_constant_10) # plot diagnostics look good

summary(mc_control_10) # model results >10% germination in control
summary(mc_dark_10) # model results >10% germination in control
summary(mc_constant_10) # model results >10% germination in control

plot(mc_control_50)  # plot diagnostics look good
plot(mc_dark_50) # plot diagnostics look good
plot(mc_constant_50) # plot diagnostics look good

summary(mc_control_50) # model results >50% germination in control
summary(mc_dark_50) # model results >50% germination in control
summary(mc_constant_50) # model results >50% germination in control

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_control_10$VCV[,"animal"]/(mc_control_10$VCV[,"animal"] + mc_control_10$VCV[,"units"]) 
mean(mc_control_10$VCV[,"animal"]/(mc_control_10$VCV[,"animal"] + mc_control_10$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_control_10$VCV[,"animal"]/(mc_control_10$VCV[,"animal"] + mc_control_10$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_control_10$VCV[,"animal"]/(mc_control_10$VCV[,"animal"] + mc_control_10$VCV[,"units"]))[, 2] %>% round(2)

lambda <- mc_dark_10$VCV[,"animal"]/(mc_dark_10$VCV[,"animal"] + mc_dark_10$VCV[,"units"]) 
mean(mc_dark_10$VCV[,"animal"]/(mc_dark_10$VCV[,"animal"] + mc_dark_10$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_dark_10$VCV[,"animal"]/(mc_dark_10$VCV[,"animal"] + mc_dark_10$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_dark_10$VCV[,"animal"]/(mc_dark_10$VCV[,"animal"] + mc_dark_10$VCV[,"units"]))[, 2] %>% round(2)

lambda <- mc_constant_10$VCV[,"animal"]/(mc_constant_10$VCV[,"animal"] + mc_constant_10$VCV[,"units"]) 
mean(mc_constant_10$VCV[,"animal"]/(mc_constant_10$VCV[,"animal"] + mc_constant_10$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_constant_10$VCV[,"animal"]/(mc_constant_10$VCV[,"animal"] + mc_constant_10$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_constant_10$VCV[,"animal"]/(mc_constant_10$VCV[,"animal"] + mc_constant_10$VCV[,"units"]))[, 2] %>% round(2)

