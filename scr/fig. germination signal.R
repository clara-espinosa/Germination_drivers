library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(lubridate);library(vegan)
library(rstatix);library(geomtextpath);library(psych)
library (ggrepel);library (ggpubr) ;library(binom)
library(RColorBrewer);library(ggtext);library(patchwork)

# germination responses
#### visualization significant differences MCMC-GLMM  fig 2 ####
x11()
# treatment
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(!(species%in%nogerm_species$species))%>%
  filter(!species%in%coldstrat_highgerm$species)%>%
  filter(!treatment=="E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>%# germinated in darkness but due to mathematical artifacts
  dplyr::select(community, habitat, species, treatment,finalgerm, viable)%>%
  group_by (treatment) %>% #, community, habitat   temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light", "Darkness"="B_alternate_dark",
                                "Constant temperatures"="D_constant_light"))%>%
  ggplot(aes(x= treatment, y= mean))+
  geom_bar(color = "black", stat = "identity", fill= "grey") +# ,position = "dodge", width = 0.9
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+ #,position= position_dodge(0.9)
  scale_y_continuous(limits = c(0,0.65))+ 
  geom_segment (aes(x= 1,xend =2,  y = 0.47, yend= 0.47), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.48, label = "·", size= 6)+
  geom_segment (aes(x= 1,xend =3,  y = 0.52, yend= 0.52), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 2, y= 0.53, label = "**", size= 6)+
  labs (x = "Treatments", y="Germination proportion", title = "Response to germination drivers", subtitle = "A) Treatment comparison (n = 38)")+
  theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.margin = margin(0,0,0,0),
         legend.key.size = unit(1,"line"),
         legend.key.spacing = unit(2, "pt"),
         axis.title.x= element_blank(), 
         axis.text.x=element_text(angle = 30, hjust = 0.9, color = "black"),
         axis.ticks.x=element_blank(),
         legend.position = "none") ->fig3a;fig3a 

# treatment x biogeography
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(!(species%in%nogerm_species$species))%>%
  filter(!species%in%coldstrat_highgerm$species)%>%
  filter(!treatment=="E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter (!species == "Cerastium ramosissimum")%>%# germinated in darkness but due to mathematical artifacts
  dplyr::select(community, habitat, biogeography, species, treatment,finalgerm, viable)%>%
  group_by (treatment, biogeography) %>% #, community, habitat   temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light", "Darkness"="B_alternate_dark",
                                "Constant temperatures"="D_constant_light"))%>%
  convert_as_factor(biogeography)%>%
  mutate(biogeography= fct_relevel(biogeography, "Mediterranean", "Eurosiberian", "Endemic", "Broad range" ))%>%
  ggplot(aes(x= treatment, y= mean, fill= biogeography))+
  geom_bar(color = "black", stat = "identity",position = "dodge", width = 0.9) +# 
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4,position= position_dodge(0.9))+
  scale_fill_manual(name = "Biogeograhical \n distribution",
                    labels = c("Mediterranean n = 10", "Eurosiberian n = 12", "Endemic n = 10", "Broad range n = 6"),
                    values = c("#FDAE61","#2B83BA","#ABDDA4","#FFFFBF"),
                    guide = guide_legend (title.position = "top",direction = "vertical")) + 
  scale_y_continuous(limits = c(0,0.65))+
  labs (x = "Treatments", y="Germination proportion", subtitle = "B) Biogeography distribution comparison")+
  theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.key.size = unit(1,"line"),
         #legend.key.spacing = unit(3, "pt"),
         legend.margin = margin(0,0,0,0),
         axis.text.x=element_text(angle = 30, hjust = 0.9, color = "black"),
         axis.title= element_blank(), 
         axis.text.y = element_blank(),
         legend.position = "right")->fig3b;fig3b




# Combine figure 
fig3a + fig3b + plot_layout(widths = c(1,1.8)) -> fig3;fig3

ggsave(filename = "germination signals.png", plot =fig3 , path = "results/Figures/", 
       device = "png", dpi = 600) #, width = 180, units = "mm"

# EXTRA ####
# community comparison
finalgerm %>%
  merge(species)%>%
  na.omit() %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>% 
  filter(!species%in%coldstrat_highgerm$species)%>%
  filter(!treatment=="E_cold_stratification")%>%
  dplyr::select(community, habitat, species, treatment,finalgerm, viable)%>%
  group_by (treatment, community) %>% #, community, habitat   temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light", "Darkness"="B_alternate_dark",
                               "Water stress"="C_alternate_WP", "Constant temperatures"="D_constant_light"))%>%
  ggplot(aes(x= treatment, y= mean, fill= treatment))+
  geom_bar(color = "black", stat = "identity") +# ,position = "dodge", width = 0.9
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+ #,position= position_dodge(0.9)
  facet_wrap(~community)+
  scale_fill_manual(name = "Treatment", values = c("#0077BC", "#41ab5d","#f6e528", "#f68b08"),
                    guide = guide_legend (title.position = "top",direction = "horizontal")) + 
  labs (x = "Treatments", y="Germination proportion", title = "Response to germination drivers")+
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 20), #hjust = 0.5,
         strip.text = element_text (size =18),
         strip.background = element_blank(), 
         panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         plot.margin =,
         legend.margin = margin(0,0,0,0),
         axis.title.x= element_blank(), 
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.position = "bottom")
