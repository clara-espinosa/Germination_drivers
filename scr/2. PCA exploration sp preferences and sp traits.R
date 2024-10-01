library (tidyverse); library (dplyr);library(rstatix); library(ggrepel)
library (ggpubr);library (vegan)
# Species levels: ecological preferences and germination responses
# Species preferences calculated from ibuttons data and inventories
# original scripts and data in Omaña80 and Picos Github repository (respectively)
#script name indices and glms

# FROM ibuttons data: in picos from 02/10/2018 to 07/08/2019; in villabandin from 12/07/2021 to 29/05/2022
    # FDD: sum of degrees (daily mean) when daily Tmean is below 0 ºC (in absolute values)
    # GDD: sum of degrees (daily mean) when daily Tmean is above 5 ºC 
    # Snowdays = sum of days with snow cover, calculated as days when Tmax <0.5 and Tmin > -0.5 ºC
# These variables (altogether with bio1, bio2 and bio7) calculated mean x day, then month and finally year (whole datset)

# From inventory data: 80 plots per each community, at each plot species coverage estimated in %
# To calculate species preferences several filters applied:
    # consider only plots where each species has at least 10% coverage
    # climatic variables weighted per coverage 


####load species data and select variables ####

read.csv("data/species.csv", sep=",") %>%
  select(species, temperature_regime, family, community, habitat)  %>%
  convert_as_factor(temperature_regime, family, community, habitat) -> species

#check species names to match; problem with minuartia CF
#setdiff(species$species, sp_pref_picos$species)

#### Species climatic preferences PCA and correlations ####
# Mediterranean
read.csv("data/sp_pref_villa.csv", sep = ",") %>%
  dplyr::select(species, community, bio1, bio2, bio7, FDD, GDD, Snw) -> sp_pref_villa

sp_pref_villa%>%
  merge(species, by = c("community", "species"))%>%
  dplyr::select(community, species, family, habitat, bio1:Snw) %>%
  mutate(Taxon = species) %>% 
  group_by (community, species, family, habitat)%>%
  summarise(across (bio1:Snw, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(species =make.cepnames(species))-> sp_pref_villa  # USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 

sp_pref_villa[, 5:10] %>%
  FactoMineR::PCA() -> pca_villa

pca_villa$var$contrib
pca_villa$eig
sp_pref_villa[, 5:10] %>%  cor()

cbind((sp_pref_villa %>%  dplyr::select(community, species, family, habitat)), data.frame(pca_villa$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds_villa

pca_villa$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars_villa

### Plot PCA
x11()
ggplot(pcaInds_villa, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars_villa, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 4, shape = 21) + # family
  #geom_label(data = pcaVars_villa, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  geom_label_repel(data = pcaVars_villa, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
                   point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds_villa, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Mediterranean species climatic preferences", tag = "A")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_villa$eig[1, 2], 0),
                                  "% variance explained)", sep = ""), limits = c(-5, 5)) + #
  scale_y_continuous(name = paste("Axis 2 (", round(pca_villa$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""), limits = c(-4, 4)) -> Med_sp_clim; Med_sp_clim
  
# Temperate 
#setdiff(sp_pref_picos$species, species$species)
#setdiff(species$species, sp_pref_picos$species)

read.csv("data/sp_pref_picos.csv", sep = ";") %>%
  select(species, community, bio1, bio2, bio7, FDD, GDD, Snw) -> sp_pref_picos

sp_pref_picos%>%
  merge(species, by = c("community", "species"))%>%
  dplyr::select(community, species, family, habitat, bio1:Snw) %>%
  mutate(Taxon = species, 
         species =make.cepnames(species))%>% # USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  group_by (community, species, family, habitat)%>%
  summarise(across (bio1:Snw, ~ mean(.x, na.rm = TRUE))) -> sp_pref_picos


sp_pref_picos[, 5:10] %>%
  FactoMineR::PCA() -> pca_picos

pca_picos$var$contrib
pca_picos$eig
sp_pref_picos[, 5:10] %>%cor()

cbind((sp_pref_picos %>%  dplyr::select(community, species, family, habitat,)), data.frame(pca_picos$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds_picos

pca_picos$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars_picos

### Plot PCA
ggplot(pcaInds_picos, aes(x = Dim.1, y = Dim.2)) +
  #coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars_picos, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 4, shape = 21) + # family
  #geom_label(data = pcaVars_picos, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  geom_label_repel(data = pcaVars_picos, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
                   point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds_picos, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Temperate species climatic preferences", tag = "B")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_picos$eig[1, 2], 0),
                                  "% variance explained)", sep = ""), limits = c(-5, 5)) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca_picos$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""), limits = c(-4, 4)) -> Tem_sp_clim; Tem_sp_clim

# combine graph
ggarrange(Med_sp_clim,Tem_sp_clim,  ncol= 2, common.legend = F)-> sp_clim;sp_clim

#### Species TRAITS PCA ###############
# Mediterranean
germ_odds_M%>% # log odds ratios
  merge(seed_mass_M)%>% # log transformed
  merge(plant_height_M)%>% # log transformed
  merge(leaves_traits_M )%>% #only leaf area log transformed
  dplyr::select(species, seed_mass,
                plant_height, 
                leaf_area, LDMC, SLA)%>% #, odds_B_dark, odds_C_WP, odds_D_constant
  merge(spe_med, by = c("species"))%>%
  dplyr::select(species, family, seed_mass:SLA) %>%
  # USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  na.omit()%>% # lose 1 species
  mutate(species =make.cepnames(species))->sp_traits_M_PCA
  


sp_traits_M_PCA[, 3:7] %>%
  FactoMineR::PCA() -> pca_traits_med

pca_traits_med$var$contrib
pca_traits_med$eig
sp_traits_M_PCA[, 3:7] %>%cor()

cbind((sp_traits_M_PCA %>%  dplyr::select(species, family)), data.frame(pca_traits_med$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds_traits_med

pca_traits_med$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable")-> pcaVars_traits_med

### Plot PCA
x11()
ggplot(pcaInds_traits_med, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars_traits_med, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 4, shape = 21) + # family
  #geom_label(data = pcaVars_picos, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  geom_label_repel(data = pcaVars_traits_med, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
                   point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds_traits_med, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Mediterranean species traits PCA", tag = "A")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_traits_med$eig[1, 2], 0),
                                  "% variance explained)", sep = ""), limits = c(-5, 5)) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca_traits_med$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""), limits = c(-4, 4)) 

# Temperate
germ_odds_T%>% # log odd ratios
  merge(seed_mass_T)%>% #log seed mass
  merge(plant_height_T)%>% # log plantheight
  merge(leaves_traits_T)%>% # only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height,
                leaf_area, LDMC, SLA)%>% #, odds_B_dark, odds_C_WP, odds_D_constant
  merge(spe_tem, by = c("species"))%>%
  dplyr::select(species, family, seed_mass:SLA) %>%
  # USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
  #na.omit()%>% # lose 3 species
  mutate(species =make.cepnames(species))->sp_traits_T_PCA

sp_traits_T_PCA[, 3:7] %>%
  FactoMineR::PCA() -> pca_traits_tem

pca_traits_tem$var$contrib
pca_traits_tem$eig
sp_traits_T_PCA[, 3:7] %>%cor()

cbind((sp_traits_T_PCA %>%  dplyr::select(species, family)), data.frame(pca_traits_tem$ind$coord[, 1:4])) %>%
  mutate(species = factor(species)) -> pcaInds_traits_tem

pca_traits_tem$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable")-> pcaVars_traits_tem

### Plot PCA
x11()
ggplot(pcaInds_traits_tem, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars_traits_tem, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = family), color = "black", show.legend = T, size = 4, shape = 21) + # family
  #geom_label(data = pcaVars_picos, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = F, size = 4) +
  geom_label_repel(data = pcaVars_traits_tem, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 5, segment.size= 1,
                   point.padding = 0.2, nudge_x = .15, nudge_y = .5,segment.curvature = -1e-20, segment.linetype = 1, segment.color = "red", arrow = arrow(length = unit(0.015, "npc")))+
  geom_text_repel (data = pcaInds_traits_tem, aes (x = Dim.1, y = Dim.2, label = species), show.legend = F, size =4) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(title= "Temperate species traits PCA", tag = "A")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.tag.position = c(0.02,1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca_traits_tem$eig[1, 2], 0),
                                  "% variance explained)", sep = ""), limits = c(-5, 5)) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca_traits_tem$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""), limits = c(-4, 4)) 
