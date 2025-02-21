library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(viridis);library(MCMCglmm); library(binom);library (ggpubr)
library(glmmTMB); library(DHARMa);library(scales)

# ANALYSIS FROM GERMINATION ECOLOGY PERSPECTIVE ####

# species with 0 germination across all experiment (remove from analysis)
  filter (!species == "Euphrasia salisburgensis")%>% # 0 germination across all treatments
  filter (!species == "Gentiana verna")%>% # 0 germination across all treatments
  filter (!species == "Gentianella campestris")%>% # 0 germination across all treatments
  filter (!species == "Kobresia myosuroides")%>% # 0 germination across all treatments
  filter (!species == "Salix breviserrata")%>% # 0 germination across all treatments
  filter (!species == "Sedum album")%>% # 0 germination across all treatments
  filter (!species == "Sedum atratum")%>% # 0 germination across all treatments

  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Teesdalia conferta")%>% # all seeds germinated under cold stratification
  filter (!species == "Phalacrocarpum oppositifolium")%>% # only germinated under cold stratification
  filter (!species == "Avenella flexuosa")%>% # only germinated under cold stratification
  # special case #
  filter (!species == "Cerastium ramosissimum")%>% # germinated under cold stratification but due to mathematical
  # artifacts with averages seems like some germination under darkness RECHECK!!!
### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL####
### Read tree
#both communities
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set priors for germination models (as many prior as random factors)
# change nu = 2 and alpha.V = 1000 for germination rate (originally nu = 1 and alpha.v = 500)
priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        #G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))
  
### 1. Cold stratification differences ####
cold_strat%>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(!(species%in%nogerm_species$species))%>%
  select(community, species,  treatment, petri, finalgerm, viable)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  convert_as_factor(species, community, petri)-> mcmc_cold

str(mcmc_cold) 
unique(mcmc_cold$species)# 50 species
mcmc_cold%>%
  group_by(species,community)%>%
  summarise(finalgerm= sum(finalgerm),
            viable=sum(viable))%>%
  group_by(community)%>%
  tally() 

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ community, # ,*community 
                   random = ~ animal + ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_cold,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_cold

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(mc_cold) # models with  both communities treatment * community (need to be updated is we want to use it)

# load("results/mcmc.Rdata")
summary(mc_cold) # need to be updated is we want to use it


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


### 2. Germination signals (Control/Light/Alternating) ####
# 2.1 average responses to germination drivers ####
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(!treatment== "E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  select(community, species, habitat, treatment, opt_temp,petri,finalgerm, viable)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>% # germinated under cold stratification but due to mathematical
  # artifacts with averages seems like some germination under darkness RECHECK!!!
  #filter(community =="Temperate")%>% #Temperate    Mediterranean
  convert_as_factor(species, community, petri,opt_temp, community)%>% 
  droplevels()-> mcmc


str(mcmc) 
unique(mcmc$species)# 38 species
unique(mcmc$treatment)

mcmc%>%
  group_by(species,community)%>%
  summarise(finalgerm= sum(finalgerm),
            viable=sum(viable))%>%
  group_by(community)%>%
  tally() 

#descriptive
mcmc%>%
  group_by(treatment)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        #G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, #
                   random = ~ animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_treat

#save(m1, file = "results/mcmc.Rdata")

x11()
plot(mc_treat) # model with all data ~treatment 
# load("results/mcmc.Rdata")
summary(mc_treat) #data ~treatment

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_treat$VCV[,"animal"]/(mc_treat$VCV[,"animal"] + mc_treat$VCV[,"units"]) 

mean(mc_treat$VCV[,"animal"]/(mc_treat$VCV[,"animal"] + mc_treat$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_treat$VCV[,"animal"]/(mc_treat$VCV[,"animal"] + mc_treat$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_treat$VCV[,"animal"]/(mc_treat$VCV[,"animal"] + mc_treat$VCV[,"units"]))[, 2] %>% round(2)


# Random effects animal
summary(mc_treat)$Gcovariances[1, 1] %>% round(2) 
summary(mc_treat)$Gcovariances[1, 2] %>% round(2) 
summary(mc_treat)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_treat)$Gcovariances[2, 1] %>% round(2)
summary(mc_treat)$Gcovariances[2, 2] %>% round(2) 
summary(mc_treat)$Gcovariances[2, 3] %>% round(2) 


# 2.2 Comparison of community within each germination treatment ####
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(treatment== "B_alternate_dark" )%>% # "B_alternate_dark"  # "D_constant_light"  "A_alternate_light"
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!species%in%coldstrat_highgerm$species)%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>% # germinated under cold stratification but due to mathematical
  # artifacts with averages seems like some germination under darkness RECHECK!!!
  convert_as_factor(species, community, petri,opt_temp, community)%>% 
  droplevels()-> mcmc_dark

str(mcmc_control) 
unique(mcmc_constant$species)# 38 species
unique(mcmc_dark$species)

mcmc_constant%>%
  group_by(species,community)%>%
  summarise(finalgerm= sum(finalgerm),
            viable=sum(viable))%>%
  group_by(community)%>%
  tally() 

#descriptive
mcmc_constant%>%
  group_by(community)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ community, #
                   random = ~ animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_constant,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_constant

#save(m1, file = "results/mcmc.Rdata")

x11()
plot(mc_control) # germination in control ~ community
plot(mc_dark) # germination in control ~ community
plot(mc_constant) # germination in control ~ community
# load("results/mcmc.Rdata")
summary(mc_control) 
summary(mc_dark)
summary(mc_constant)

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_dark$VCV[,"animal"]/(mc_dark$VCV[,"animal"] + mc_dark$VCV[,"units"]) 

mean(mc_dark$VCV[,"animal"]/(mc_dark$VCV[,"animal"] + mc_dark$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_dark$VCV[,"animal"]/(mc_dark$VCV[,"animal"] + mc_dark$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_dark$VCV[,"animal"]/(mc_dark$VCV[,"animal"] + mc_dark$VCV[,"units"]))[, 2] %>% round(2)


# Random effects animal
summary(mc_treat)$Gcovariances[1, 1] %>% round(2) 
summary(mc_treat)$Gcovariances[1, 2] %>% round(2) 
summary(mc_treat)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_treat)$Gcovariances[2, 1] %>% round(2)
summary(mc_treat)$Gcovariances[2, 2] %>% round(2) 
summary(mc_treat)$Gcovariances[2, 3] %>% round(2) 

### 3. Water stress differences ####
# Response to water stress calculated as the percentual reduction compared to control.
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>% # germinated under cold stratification but due to mathematical artifacts
  mutate(germpro = finalgerm/viable)%>%
  group_by(species)%>%
  summarise(germpro=mean(germpro))%>%
  filter(germpro==0)-> nogerm_control_waterstress

finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>% # germinated under cold stratification but due to mathematical artifacts
  filter(!(species%in%nogerm_control_waterstress$species))%>%
  convert_as_factor(species, community, petri, community)%>% 
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  droplevels()-> mcmc_water

finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>% # germinated under cold stratification but due to mathematical artifacts
  filter(!(species%in%nogerm_control_waterstress$species))%>%
  mutate(germpro = finalgerm/viable)%>% # generate Nan because of 0 viable
  mutate_all(~replace(., is.nan(.), 0))%>%
  select(community, species, habitat, community, treatment, petri,germpro)%>%
  spread(treatment, germpro)%>%
  mutate(germ_reduction=(1-(C_alternate_WP/A_alternate_light))*100)%>%
  mutate_all(~replace(., is.nan(.), 0))%>% # generate Nan because 0 germ pro
  convert_as_factor(species, community, petri, community)%>% 
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  #mutate(germ_reduction=sqrt(germ_reduction))%>% problem with negative values in those cases were germination in water stress was slightly higher
  droplevels()-> mcmc_waterreduction

str(mcmc_water) 
unique(mcmc_water$species)# 38 species

str(mcmc_waterreduction) 
unique(mcmc_waterreduction$species)

mcmc_waterreduction%>%
  group_by(species,community)%>%
  summarise(germ_reduction= mean(germ_reduction))%>%
  print(n= 34)
  group_by(community)%>%
  tally() 

#descriptive
mcmc_water_reduction%>%
  group_by(species, community)%>%
  summarise(germ_reduction= mean(germ_reduction))%>%
  group_by(community)%>%
  get_summary_stats (germ_reduction)

#descriptive
mcmc_water%>%
  group_by(treatment)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, #
                   random = ~ animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_water,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_water


MCMCglmm::MCMCglmm(germ_reduction ~ community, #
                   random = ~ animal + ID , #
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = mcmc_waterreduction,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_water_reduction
#save(m1, file = "results/mcmc.Rdata")

x11()
plot(mc_water) # germination ~ treatment
plot(mc_water_reduction) # germination reduction ~community
# load("results/mcmc.Rdata")
summary(mc_water) 
summary(mc_water_reduction) 

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"]) 

mean(mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"]))[, 2] %>% round(2)


# Random effects animal
summary(mc_treat)$Gcovariances[1, 1] %>% round(2) 
summary(mc_treat)$Gcovariances[1, 2] %>% round(2) 
summary(mc_treat)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_treat)$Gcovariances[2, 1] %>% round(2)
summary(mc_treat)$Gcovariances[2, 2] %>% round(2) 
summary(mc_treat)$Gcovariances[2, 3] %>% round(2) 