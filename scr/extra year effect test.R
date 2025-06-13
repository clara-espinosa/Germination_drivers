library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(binom)


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
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### 1. Cold stratification differences ####
cold_strat%>%
  merge(read.csv("data/species.csv", sep = ";"), by = c("species", "code"))%>%
  filter(!(species%in%nogerm_species$species))%>%
  select(community, species,  treatment, petri, finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  convert_as_factor(species, community, petri, collection_year)-> mcmc_cold

str(mcmc_cold) 
unique(mcmc_cold$species)# 50 species
mcmc_cold%>%
  group_by(collection_year,community)%>%
  summarise(finalgerm= sum(finalgerm),
            viable=sum(viable))%>%
  group_by(community)%>%
  tally() 

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ community, # ,*community 
                   random = ~ collection_year + animal + ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_cold,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_cold2

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(mc_cold2) # collection_year included as random
# load("results/mcmc.Rdata")
summary(mc_cold2) # collection_year included as random

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_cold2$VCV[,"animal"]/(mc_cold2$VCV[,"animal"] + mc_cold2$VCV[,"units"]) 

mean(mc_cold2$VCV[,"animal"]/(mc_cold2$VCV[,"animal"] + mc_cold2$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_cold2$VCV[,"animal"]/(mc_cold2$VCV[,"animal"] + mc_cold2$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_cold2$VCV[,"animal"]/(mc_cold2$VCV[,"animal"] + mc_cold2$VCV[,"units"]))[, 2] %>% round(2)

### 2. Fine scale germination signals (Control/Light/Alternating) ####
# 2.1 average responses to germination cues compared to control ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep = ";"), by = c("species", "code"))%>%
  filter(!treatment== "E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  select(community, species, habitat, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  #filter(community =="Temperate")%>% #Temperate    Mediterranean
  convert_as_factor(species, community, petri,opt_temp, community, collection_year)%>% 
  droplevels()-> mcmc

#descriptive
mcmc%>%
  group_by(treatment)%>%
  summarise(finalgerm= sum(finalgerm), viable= (sum(viable)))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, #
                   random = ~ collection_year+animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_treat4

#save(m1, file = "results/mcmc.Rdata")

x11()
plot(mc_treat2) # with optimal temperatures as random
plot(mc_treat3)# with optimal temperatures as fixed
plot(mc_treat4)# with collection year included as random
# load("results/mcmc.Rdata")
summary(mc_treat2) # optimal temp not significantly explains variation
summary(mc_treat3) # optimal temp not affect germination 
summary(mc_treat4) # with collection year included as random

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_treat4$VCV[,"animal"]/(mc_treat4$VCV[,"animal"] + mc_treat4$VCV[,"units"]) 

mean(mc_treat4$VCV[,"animal"]/(mc_treat4$VCV[,"animal"] + mc_treat4$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_treat4$VCV[,"animal"]/(mc_treat4$VCV[,"animal"] + mc_treat4$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_treat4$VCV[,"animal"]/(mc_treat4$VCV[,"animal"] + mc_treat4$VCV[,"units"]))[, 2] %>% round(2)


# 2.2 Comparison of community within each germination treatment ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep= ";"), by = c("species", "code"))%>%
  filter(treatment== "A_alternate_light" )%>% # "B_alternate_dark"  # "D_constant_light"  "A_alternate_light"
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!species%in%coldstrat_highgerm$species)%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  convert_as_factor(species, community, petri,opt_temp, community, collection_year)%>% 
  droplevels()-> mcmc_control
str(mcmc_control) 
unique(mcmc_constant$species)# 39 species
unique(mcmc_dark$species)


### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ community, #
                   random = ~ collection_year+animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_dark,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_dark2

#save(m1, file = "results/mcmc.Rdata")

x11()
plot(mc_control2) # germination in control ~ community, collection year included as random
plot(mc_dark2) # germination in dark ~ community, collection year included as random
plot(mc_constant2) # germination constant temp ~ community, collection year included as random
# load("results/mcmc.Rdata")
summary(mc_control2) # germination in control ~ community, collection year included as random
summary(mc_dark2) # germination in dark ~ community, collection year included as random
summary(mc_constant2) # germination constant temp ~ community, collection year included as random

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_constant$VCV[,"animal"]/(mc_constant$VCV[,"animal"] + mc_constant$VCV[,"units"]) 

mean(mc_constant$VCV[,"animal"]/(mc_constant$VCV[,"animal"] + mc_constant$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_constant$VCV[,"animal"]/(mc_constant$VCV[,"animal"] + mc_constant$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_constant$VCV[,"animal"]/(mc_constant$VCV[,"animal"] + mc_constant$VCV[,"units"]))[, 2] %>% round(2)


### 3. Water limitation differences ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep = ";"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  mutate(germpro = finalgerm/viable)%>%
  group_by(species)%>%
  summarise(germpro=mean(germpro))%>%
  filter(germpro==0)-> nogerm_control_waterstress # species with no germ either in control or water stress, remove from further analysis


#3.1 Control vs water limitation germination ####
finalgerm %>%
  merge(read.csv("data/species.csv", sep = ";"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter(!(species%in%nogerm_control_waterstress$species))%>%
  convert_as_factor(species, community, petri, community, collection_year)%>% 
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  droplevels()-> mcmc_water
str(mcmc_water) 
unique(mcmc_water$species)# 33 species


### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, #
                   random = ~ collection_year + animal + ID , #
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_water,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_water2

x11()
plot(mc_water2) # germination ~ treatment, collection year included as random
summary(mc_water2) # collection year included as random
### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_water$VCV[,"animal"]/(mc_water$VCV[,"animal"] + mc_water$VCV[,"units"]) 

mean(mc_water$VCV[,"animal"]/(mc_water$VCV[,"animal"] + mc_water$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_water$VCV[,"animal"]/(mc_water$VCV[,"animal"] + mc_water$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_water$VCV[,"animal"]/(mc_water$VCV[,"animal"] + mc_water$VCV[,"units"]))[, 2] %>% round(2)

#3.2 germiantion reduction (%) due to water stress compared to control ####
# Response to water stress calculated as the percentual reduction compared to control.
finalgerm %>%
  merge(read.csv("data/species.csv", sep=";"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable, collection_year)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter(!(species%in%nogerm_control_waterstress$species))%>%
  mutate(germpro = finalgerm/viable)%>% # generate Nan because of 0 viable
  mutate_all(~replace(., is.nan(.), 0))%>%
  select(community, species, habitat, community, treatment, petri,germpro, collection_year)%>%
  spread(treatment, germpro)%>%
  mutate(germ_reduction=(1-(C_alternate_WP/A_alternate_light))*100)%>%
  mutate_all(~replace(., is.nan(.), 0))%>% # generate Nan because 0 germ pro
  convert_as_factor(species, community, petri, community, collection_year)%>% 
  #mutate(germ_reduction=ifelse(germ_reduction<0, 0, germ_reduction))%>% #problem with negative values in those cases were germination in water stress was slightly higher
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  droplevels()%>% 
  as.data.frame()-> mcmc_waterreduction2

str(mcmc_waterreduction) # with gaussian family
str(mcmc_waterreduction2) # with exponential family
unique(mcmc_waterreduction$species)
hist(mcmc_waterreduction$germ_reduction)

### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
#)))

MCMCglmm::MCMCglmm(germ_reduction ~ community, #
                   random = ~ collection_year + animal + ID , #
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = mcmc_waterreduction2,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_water_reduction4

plot(mc_water_reduction3) # collection_year included as random Gaussian family
plot(mc_water_reduction4)# germination reduction ~community EXPONENTIAL family # collection_year included as random

summary(mc_water_reduction3)# collection_year included as randomGaussian family
summary(mc_water_reduction4) # germination reduction ~community EXPONENTIAL family # collection_year included as random

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"]) 

mean(mc_water_reduction2$VCV[,"animal"]/(mc_water_reduction2$VCV[,"animal"] + mc_water_reduction$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_water_reduction$VCV[,"animal"]/(mc_water_reduction$VCV[,"animal"] + mc_water_reduction$VCV[,"units"]))[, 2] %>% round(2)
