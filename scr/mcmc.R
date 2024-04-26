library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(viridis);library(MCMCglmm); library(binom)
library(glmmTMB); library(DHARMa);library(scales)

# ANALYSIS FROM GERMINATION ECOLOGY PERSPECTIVE ####
# dataframe with raw data ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) -> raw_df
str(raw_df)
unique(raw_df$species)

# species with 0 germination across all experiment (remove from analysis)
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")

# dataframe with species data ####
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat)  %>%
  convert_as_factor(species, code, family, community, habitat)-> species 

# viables x petri x treatment ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  select (species, code, treatment, petri, viable)%>%
  arrange(species, code, treatment, petri)-> viables_petri

# viables x species x treatment
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  select (species, code, treatment, petri, viable)%>%
  group_by (species, code, treatment)%>%
  summarise(viable= sum(viable))%>%
  arrange(species, treatment) -> viable_sp

# D0 germination check change to cold_stratification treatment and mean values ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, petri, initial, cold_strat, empty, fungus) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  mutate (viable = initial -(empty + fungus)) %>%
  group_by (species, code, petri) %>%
  summarise (viable = round(mean(viable),0),
             finalgerm= round(mean(cold_strat),0)) %>%
  mutate (finalgerm = replace_na(finalgerm , 0)) %>% # salix and solidago have some petridish with all fungus
  mutate (treatment = "E_cold_stratification")%>%
  group_by (species, code, treatment, petri) %>%
  summarize(finalgerm = round(mean (finalgerm),0), 
            viable = round(mean(viable),0))-> cold_strat

#final germination percentage calculation ####
read.csv("data/raw_data.csv", sep = ",") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species, code, treatment, petri)%>%
  summarize (finalgerm = sum(germ)) %>% # 
  select(species, code, treatment, petri, finalgerm)%>% 
  merge(viables_petri) %>%
  rbind(cold_strat)%>%
  arrange(species, code, treatment, petri)%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code, treatment, petri, finalgerm, viable,  mean, upper, lower)-> finalgerm # final germ per petri dish 

# species preferences (focus on GDD, FDD) ####
# Mediterranean
read.csv("data/sp_pref_villa.csv", sep = ",") %>%
  dplyr::select(species, community, bio1, bio2, bio7, FDD, GDD, Snw) -> sp_pref_villa
# Temperate 
read.csv("data/sp_pref_picos.csv", sep = ";") %>%
  select(species, community, bio1, bio2, bio7, FDD, GDD, Snw) -> sp_pref_picos

sp_pref_villa%>%
  rbind(sp_pref_picos)->sp_pref

species%>%
  merge(sp_pref, by =c("species", "community"))%>%
  select(species, community, code, habitat, GDD, FDD)->species2


# MCMC GLMM ####
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>%
  merge(species2)%>%
  select(community, species, habitat, FDD, GDD, treatment, finalgerm, viable)%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  na.omit ()%>%
  #filter(community =="Temperate")%>% #Temperate    Mediterranean
  as.data.frame()-> mcmc_ibc
str(mcmc_ibc)
unique(mcmc_ibc$habitat)
#### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL####
### Read tree
phangorn::nnls.tree(cophenetic(ape::read.tree("results/treegerm.tree")), 
                    ape::read.tree("results/treegerm.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

# shorter iterations
# nite = 10000
#nthi = 10
#nbur = 100

### Set priors for germination models (as many prior as random factors)
# change nu = 2 and alpha.V = 1000 for germination rate (originally nu = 1 and alpha.v = 500)
priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))
#G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))   
### TEST
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment*GDD + treatment*FDD, # *community, habitat, GDD, FDD
                   random = ~ animal + ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_ibc,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)


# load("results/mcmc.Rdata")
summary(m1)

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]) 

mean(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"])) %>% round(2)
coda::HPDinterval(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(m1)$Gcovariances[1, 1] %>% round(2) 
summary(m1)$Gcovariances[1, 2] %>% round(2) 
summary(m1)$Gcovariances[1, 3] %>% round(2)

# Random effects code:ID
summary(m1)$Gcovariances[2, 1] %>% round(2)
summary(m1)$Gcovariances[2, 2] %>% round(2) 
summary(m1)$Gcovariances[2, 3] %>% round(2) 


#### visualization significant differences MCMC-GLMM ####
# treatment x community, only differnce in WP, lower in temperate
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>% #-> mcmc_ibc
  merge(species2)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  na.omit() %>%
  group_by (treatment, community) %>% #, temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  #filter(!treatment =="B_alternate_dark")%>%
  #filter(!treatment =="C_alternate_WP")%>%
  #filter(!treatment =="D_constant_light")%>%
  ggplot()+
  geom_bar(aes(x= treatment, y= mean, fill= treatment), color = "black", stat = "identity") +
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", size=1, width = 0.4)+
  facet_grid(~community) +
  ylim (0,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "dimgrey", "chocolate1","#7AD151FF")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "Treatments", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =20),
         axis.title.x= element_blank(), 
         axis.text.x= element_blank(), 
         legend.position = "bottom")
show_col(viridis(4))

# treatment x habitat
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>% #-> mcmc_ibc
  merge(species2)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  na.omit() %>%
  group_by (treatment, habitat) %>% #, temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  #filter(!treatment =="B_alternate_dark")%>%
  #filter(!treatment =="C_alternate_WP")%>%
  #filter(!treatment =="D_constant_light")%>%
  ggplot()+
  geom_bar(aes(x= treatment, y= mean, fill= treatment), color = "black", stat = "identity") +
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", size=1, width = 0.4)+
  facet_grid(~habitat) +
  ylim (0,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "dimgrey", "chocolate1","#7AD151FF")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "Treatments", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =20),
         axis.title.x= element_blank(), 
         axis.text.x= element_blank(), 
         legend.position = "bottom")

# GDD/FDD correlation with treatment germination responses
 ##GDD        
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>% #-> mcmc_ibc
  merge(species2)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  ggplot(aes(x= GDD, y= mean, fill= treatment), color = "black")+
  geom_point(stat = "identity", shape=21, size=4) +
  facet_grid(~treatment) + 
  geom_smooth(method = "loess")+
  ylim (-0.1,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "dimgrey", "chocolate1","#7AD151FF")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "GDD (Growing dregree days)", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =16),
         axis.title.x= element_text (size =18), 
         axis.text.x= element_blank(), 
         legend.position = "bottom")
##FDD
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>% #-> mcmc_ibc
  merge(species2)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  ggplot(aes(x= FDD, y= mean, fill= treatment), color = "black")+
  geom_point(stat = "identity", shape=21, size=4) +
  facet_grid(~treatment) +
  geom_smooth(method = "lm")+
  ylim (-0.1,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "dimgrey", "chocolate1","#7AD151FF")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "FDD (Freezing degree days)", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =16),
         axis.title.x= element_text (size =18), 
         axis.text.x= element_blank(), 
         legend.position = "none")

