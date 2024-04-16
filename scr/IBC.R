# IBC exploration
library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(MCMCglmm); library(binom)
library(glmmTMB); library(DHARMa);library(scales)

# dataframe with raw data ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) -> raw_df
str(raw_df)
unique(raw_df$species)

# dataframe with species data ####
read.csv("data/species.csv", sep=";") %>%
  select(species, code, temperature_regime, family, community, habitat)  %>%
  convert_as_factor(species, code, temperature_regime, family, community, habitat) %>% #-> species
  select(species, family, community)%>%
  distinct()-> species2

# species with 0 germination across all experiment (remove from analysis)
filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")

# viables x community x treatment
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  mutate (viablePER = (viable/initial)*100) %>%
  select (species, code, treatment, petri, viable, viablePER)%>%
  merge(species)%>%
  group_by (community, treatment)%>%
  summarise(viable= sum(viable))-> viables_community
# viables x species x treatment
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  mutate (viablePER = (viable/initial)*100) %>%
  select (species, code, treatment, petri, viable)%>%
  rbind (cold_strat_viables)%>% #->viables_petri
  arrange(species, treatment)%>%
  merge(species)%>%
  group_by (community, species,treatment)%>%
  summarise(viable= sum(viable))%>%
  arrange(species, treatment) -> viable_sp

# D0 germination check change to cold_stratification treatment and mean values ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  merge(viables_petri) %>%
  select (species, code, petri, viable, initial, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code,petri) %>%
  summarize(finalgerm = sum(cold_strat),
            viable = sum(viable)) %>% 
  #mutate (viablePER = (viable/initial)*100) %>%
  mutate (finalgerm = replace_na(finalgerm , 0)) %>% # salix and solidago have some petridish with all fungus
  mutate (treatment = "E_cold_stratification") %>%
  #select(species, code, treatment, petri, viable)->cold_strat_viables
  mutate(treatment = as.factor(treatment)) %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code,  treatment, finalgerm, viable, mean, upper, lower)-> cold_strat

cold_strat%>%
  merge (species)%>%
  group_by (community, species, treatment) %>%
  summarize(finalgerm = sum (finalgerm))-> cold_strat2

# final germination
#final germination percentage calculation discounting cold stratification germ ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species, code, treatment, petri)%>%
  summarize (finalgerm = sum(germ), 
             initial = first(initial)) %>% # 
  merge(viables_petri) %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code, treatment, petri, finalgerm, viable,  mean, upper, lower)-> finalgerm # final germ per petri dish 

#visualization bars + error bars x community/treatment
finalgerm %>%
  merge(species)%>%
  group_by (community, species, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  rbind(cold_strat2)%>%
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  merge(viable_sp, by= c("community", "species", "treatment"))%>% #-> mcmc_ibc
  arrange(species, treatment) %>%
  #merge (species, by = "species") %>%
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
  filter(!treatment =="E_cold_stratification")%>%
  filter(!treatment =="B_alternate_dark")%>%
  filter(!treatment =="C_alternate_WP")%>%
  #filter(!treatment =="D_constant_light")%>%
  ggplot()+
  geom_bar(aes(x= treatment, y= mean, fill= treatment), color = "black", stat = "identity") +
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", size=1, width = 0.4)+
  facet_grid(~community) +
  ylim (0,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Constant Temperature"), values = c("#2A788EFF", "#7AD151FF")) + 
  #"Darkness" = dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "Treatments", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =20),
         axis.title.x= element_blank(), 
         axis.text.x= element_blank(), 
         legend.position = "none")
show_col(viridis(4))
#MCMC GLMM IBC
finalgerm %>%
  merge(species)%>%
  group_by (community, species, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  rbind(cold_strat2)%>%
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  merge(viable_sp, by= c("community", "species", "treatment"))%>%
  merge(species2)%>%
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
  select(!family) %>%
  filter(community =="Temperate")%>% #Temperate    Mediterranean
  as.data.frame()-> mcmc_ibc
unique(mcmc_ibc$species)

a <- glmmTMB(cbind(finalgerm, viable - finalgerm) ~ treatment*community,  family = binomial, data= mcmc_ibc) 
summary(a)
#separate communities
a <- glmmTMB(cbind(finalgerm, viable - finalgerm) ~ treatment,  family = binomial, data= mcmc_ibc) 
summary(a)

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
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, # *community
                   random = ~ animal + ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_ibc,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)


# load("results/mcmc.Rdata")
summary(m1)

