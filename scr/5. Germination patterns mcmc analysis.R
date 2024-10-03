library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(viridis);library(MCMCglmm); library(binom);library (ggpubr)
library(glmmTMB); library(DHARMa);library(scales)

# ANALYSIS FROM GERMINATION ECOLOGY PERSPECTIVE (as germination phenology) ####

# species with 0 germination across all experiment (remove from analysis)
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")%>%
  filter (!species == "Teesdalia conferta")%>%

# dataframe with species data ####
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat)  %>%
  convert_as_factor(species, code, family, community, habitat)-> species 
 unique(species$species)  

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
  arrange(species, treatment) -> viables_sp

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
x11()
cold_strat%>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm), viable = sum(viable)) %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  merge (species)%>%
  group_by (community)%>%
  summarize(finalgerm = sum (finalgerm), viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  ggplot()+
  geom_bar(aes(x= community, y= mean), fill= "dimgrey", color = "black", stat = "identity") +
  geom_errorbar(aes(x= community, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+
  ylim (0,0.5)+
  labs (x = "Community", y="Germination proportion", title = "Cold stratification", 
        subtitle = "Pretreatment", tag = "A")+ #
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 20), #hjust = 0.5,
         axis.text.x = element_text (color="black"),
         legend.position = "bottom") -> fig2a;fig2a

#final germination percentage calculation ####
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
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))-> finalgerm 

# MCMC GLMM ####
finalgerm %>%
  merge(species)%>%
  select(community, species, habitat, treatment, finalgerm, viable)%>%
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
  filter (!species == "Teesdalia conferta")%>%
  filter(community =="Temperate")%>% #Temperate    Mediterranean
  na.omit ()-> mcmc_tem


str(mcmc)# all species
str (mcmc_tem) # temperate species
str(mcmc_med) # mediterranean species
unique(mcmc$species)
unique(mcmc_med$species)
unique(mcmc_tem$species)
#### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL####
### Read tree
#both communities
phangorn::nnls.tree(cophenetic(ape::read.tree("results/treegerm.tree")), 
                    ape::read.tree("results/treegerm.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
# ONly temperate community
phangorn::nnls.tree(cophenetic(ape::read.tree("results/TEM_tree.tree")), 
                    ape::read.tree("results/TEM_tree.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
#ONLY mediterranean community
phangorn::nnls.tree(cophenetic(ape::read.tree("results/MED_tree.tree")), 
                    ape::read.tree("results/MED_tree.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL
 plot(ape::read.tree("results/MED_tree.tree"))
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
MCMCglmm::MCMCglmm(cbind(finalgerm, viable - finalgerm) ~ treatment, # ,*community 
                   random = ~ animal + ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = mcmc_med,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> mc_med

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(mc_both) # models with  both communities treatment * community
plot(mc_tem) 
plot(mc_med)

# load("results/mcmc.Rdata")
summary(mc_both)
summary(mc_tem)
summary(mc_med)

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- mc_med$VCV[,"animal"]/(mc_med$VCV[,"animal"] + mc_med$VCV[,"units"]) 

mean(mc_med$VCV[,"animal"]/(mc_med$VCV[,"animal"] + mc_med$VCV[,"units"])) %>% round(2)
coda::HPDinterval(mc_med$VCV[,"animal"]/(mc_med$VCV[,"animal"] + mc_med$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(mc_med$VCV[,"animal"]/(mc_med$VCV[,"animal"] + mc_med$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(mc_med)$Gcovariances[1, 1] %>% round(2) 
summary(mc_med)$Gcovariances[1, 2] %>% round(2) 
summary(mc_med)$Gcovariances[1, 3] %>% round(2)

# Random effects ID
summary(mc_med)$Gcovariances[2, 1] %>% round(2)
summary(mc_med)$Gcovariances[2, 2] %>% round(2) 
summary(mc_med)$Gcovariances[2, 3] %>% round(2) 


#### visualization significant differences MCMC-GLMM ####
# treatment x community, only differnce in WP, lower in temperate
ann_text_m <- data.frame (x=1.6, y= 0.41, label = "*", community = "Mediterranean")
ann_text_t <- data.frame (x=1.6, y= 0.41, label = "*", community = "Temperate")
x11()
finalgerm %>%
  merge(species)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  filter(!species == "Teesdalia conferta") %>%
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
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+
  facet_grid(~community) +
  ylim (0,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "#440154FF", "#FDE725FF","#35B779FF"),
                    guide = guide_legend (title.position = "top",direction = "horizontal")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "Treatments", y="Germination proportion", title = "Response to germination drivers", tag = "B")+ #
  # add significances
  geom_segment (aes(x= 1,xend =3.1,  y = 0.45, yend= 0.45), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 2, y= 0.46, label = "***", size= 6)+
  geom_segment(aes(x=1, xend =2.1, y= 0.4, yend=0.4), color = "black", linewidth = 1.3,inherit.aes=FALSE)+
  geom_text(data=ann_text_m, label ="***", aes(x=x, y=y), size= 6)+
  geom_text(data=ann_text_t, label ="*", aes(x=x, y=y), size= 6)+
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 20), #hjust = 0.5,
         strip.text = element_text (size =18),
         strip.background = element_blank(), 
         panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         plot.margin =,
         legend.margin = margin(0,0,0,0),
         axis.title.x= element_blank(), 
         axis.text.x= element_blank(), 
         legend.position = c(0.5,-0.08))-> fig2b;fig2b
show_col(viridis(4))

# Combine figure 2 
library(patchwork)
fig2a + fig2b + plot_layout(widths = c(1,2))->Fig2;Fig2 

ggsave(filename = "Fig2.png", plot =Fig2 , path = "results/Figures/", 
       device = "png", dpi = 600)
