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
  dplyr::select(community, habitat, species, treatment,finalgerm, viable)%>%
  group_by (treatment) %>% #, community, habitat   temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light", "Darkness"="B_alternate_dark",
                                "Constant \n temperatures"="D_constant_light"))%>%
  mutate(treatment = fct_relevel(treatment, "Control", "Darkness", "Constant \n temperatures"))%>%
  ggplot(aes(x= treatment, y= mean))+
  geom_bar(color = "black", stat = "identity", fill= "grey") +# ,position = "dodge", width = 0.9
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+ #,position= position_dodge(0.9)
  scale_y_continuous(limits = c(0,0.6))+ 
  geom_segment (aes(x= 1,xend =2,  y = 0.47, yend= 0.47), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.48, label = "***", size= 6)+
  geom_segment (aes(x= 1,xend =3,  y = 0.52, yend= 0.52), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 2, y= 0.53, label = "*", size= 6)+
  labs ( x = "Treatments", y="Germination proportion",  subtitle = "A) By Treatment")+
  theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.margin = margin(0,0,0,0),
         legend.key.size = unit(1,"line"),
         legend.key.spacing = unit(2, "pt"),
         axis.title.x= element_blank(), 
         axis.text.x=element_text(color = "black"), #angle = 30, hjust = 0.9, 
         axis.ticks.x=element_blank(),
         legend.position = "none") ->fig4a;fig4a 

# treatment x community
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(!(species%in%nogerm_species$species))%>%
  filter(!(species%in%coldstrat_highgerm$species))%>%
  filter(!treatment=="E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  dplyr::select(community, habitat, community, species, treatment,finalgerm, viable)%>%
  group_by (treatment, community) %>% #, community, habitat   temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light", "Darkness"="B_alternate_dark",
                                "Constant \n temperatures"="D_constant_light"))%>%
  convert_as_factor(community)%>%
  mutate(community= fct_relevel(community, "Mediterranean", "Temperate"))%>%
  ggplot(aes(x= treatment, y= mean, fill= community))+
  geom_bar(color = "black", stat = "identity",position = "dodge", width = 0.9) +# 
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4,position= position_dodge(0.9))+
  scale_fill_manual(name = "Habitat",
                    labels = c("Mediterranean \n Alpine (n = 16)", "Temperate \n Alpine (n = 24)"),
                    values = c("darkgoldenrod1", "forestgreen"),
                    guide = guide_legend (title.position = "top",direction = "vertical")) + 
  scale_y_continuous(limits = c(0,0.6))+
  labs (x = "Treatments", y="Germination proportion", subtitle = "B) By Habitat")+
  theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.key.size = unit(1,"line"),
         #legend.key.spacing = unit(5, "pt"),
         legend.margin = margin(0,0,0,0),
         axis.text.x=element_text( color = "black"), #angle = 30, hjust = 0.9,
         axis.title= element_blank(), 
         axis.text.y = element_blank(),
         legend.position = "right")->fig4b;fig4b



# Combine figure 
library(patchwork)
fig4a + fig4b + plot_layout(widths = c(1.2,1.8)) + plot_annotation (title = "Response to dark and constant temperatures (n = 39)") -> fig4;fig4

ggsave(filename = "germination cues.png", plot =fig4 , path = "results/Figures/", 
       device = "png", dpi = 600,height = 150, width = 170, units = "mm") #,
