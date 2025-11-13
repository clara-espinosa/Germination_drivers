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
  #filter(!(species%in%nogerm_species$species))%>%
  #filter(!species%in%coldstrat_highgerm$species)%>%
  filter(!treatment=="E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  #filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter(species%in%above10_control$species)%>%
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
  geom_bar(color = "black", stat = "identity", fill= "grey", linewidth=0.3) +# ,position = "dodge", width = 0.9
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=0.3, width = 0.5)+ #,position= position_dodge(0.9)
  scale_y_continuous(limits = c(0,0.67))+ 
  geom_segment (aes(x= 1,xend =2,  y = 0.63, yend= 0.63), color = "black", linewidth = 0.3, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.64, label = "***", size= 2)+
  labs ( x = "Treatments", y="Germination proportion",  subtitle = "a) By Treatment")+
  theme_classic (base_size = 6) + #theme_minimal for all species for mean treatment
  theme (panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.margin = margin(0,0,0,0),
         legend.key.size = unit(0.3,"line"),
         legend.key.spacing = unit(2, "pt"),
         axis.title.x= element_blank(), 
         axis.text.x=element_text(color = "black"), #angle = 30, hjust = 0.9, 
         axis.ticks.x=element_blank(),
         legend.position = "none") ->fig4a;fig4a 

# treatment x community
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  #filter(!(species%in%nogerm_species$species))%>%
  #filter(!(species%in%coldstrat_highgerm$species))%>%
  filter(!treatment=="E_cold_stratification")%>%
  filter(!treatment== "C_alternate_WP")%>%
  #filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter(species%in%above10_control$species)%>%
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
  geom_bar(color = "black", stat = "identity",position = "dodge", width = 0.9, linewidth=0.3) +# 
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=0.3, width = 0.5,position= position_dodge(0.9))+
  scale_fill_manual(name = "Habitat",
                    labels = c("Mediterranean \n Alpine (n = 13)", "Temperate \n Alpine (n = 17)"),
                    values = c("darkgoldenrod1", "forestgreen"),
                    guide = guide_legend (title.position = "left",direction = "horizontal")) + 
  scale_y_continuous(limits = c(0,0.67))+
  labs (x = "Treatments", y="Germination proportion", subtitle = "b) By Habitat")+
  theme_classic (base_size = 6) + #theme_minimal for all species for mean treatment
  theme (panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.key.size = unit(1,"line"),
         legend.margin = margin(0,0,0,0),
         axis.text.x=element_text( color = "black"), #angle = 30, hjust = 0.9,
         axis.title= element_blank(), 
         axis.text.y = element_blank(),
         legend.position = "bottom")->fig4b;fig4b



# Combine figure 
library(patchwork)
fig4a + fig4b + plot_layout(widths = c(1.2,1.8), guides = "collect")& theme(legend.position = "bottom") -> fig4;fig4

ggsave(filename = "Fig5.png", plot =fig4 , path = "results/Figures/", 
       device = "png", dpi = 600, width = 80, height = 60, units = "mm") #,
