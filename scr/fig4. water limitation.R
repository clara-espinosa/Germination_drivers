library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(lubridate);library(vegan)
library(rstatix);library(geomtextpath);library(psych)
library (ggrepel);library (ggpubr) ;library(binom)
library(RColorBrewer);library(ggtext);library(patchwork)

# Water stress responses visualization
### A) differences between control and wp ####
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, treatment, opt_temp,petri,finalgerm, viable)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea") %>% # only germinated under cold stratification
  filter(!(species%in%nogerm_control_waterstress$species))%>%
  group_by (treatment) %>% #, community, habitat   temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(treatment = as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light", "Water limitation"="C_alternate_WP"))%>%
  ggplot(aes(x= treatment, y= mean))+
  geom_bar(color = "black", stat = "identity", fill= "grey") +# ,position = "dodge", width = 0.9
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=0.5, width = 0.4)+ #,position= position_dodge(0.9)
  scale_y_continuous(limits = c(0,0.6))+ 
  geom_segment (aes(x= 1,xend =2,  y = 0.54, yend= 0.54), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.55, label = "***", size= 4)+
  labs (x = "Treatments", y="Germination proportion",  subtitle = "a) Germination by Treatment")+
  theme_classic (base_size = 9) + #theme_minimal for all species for mean treatment
  theme (panel.background = element_blank(),
         legend.title = element_text(hjust=0.5),
         legend.margin = margin(0,0,0,0),
         legend.key.size = unit(1,"line"),
         legend.key.spacing = unit(2, "pt"),
         axis.title.x= element_blank(), 
         axis.text.x=element_text(color = "black"),#angle = 30, hjust = 0.9, 
         axis.ticks.x=element_blank(),
         legend.position = "none")->fig3a;fig3a

### B) differences between community ####
x11()
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter(!(species%in%nogerm_control_waterstress$species)) %>%
  mutate(germpro = finalgerm/viable)%>% # generate Nan because of 0 viable
  mutate_all(~replace(., is.nan(.), 0))%>%
  select(community, species, habitat, community, treatment, petri,germpro)%>%
  spread(treatment, germpro)%>%
  mutate(germ_reduction=(1-(C_alternate_WP/A_alternate_light))*100)%>%
  mutate_all(~replace(., is.nan(.), 0))%>% # generate Nan because 0 germ pro
  convert_as_factor(species, community, petri, community)%>%
  mutate(community= fct_relevel(community, "Mediterranean", "Temperate" ))%>%
  mutate(community= fct_recode(community,"Mediterranean \n Alpine n=15" ="Mediterranean", 
                               "Temperate \n Alpine n=19"="Temperate" ))%>%
  group_by(community)%>%
  get_summary_stats (germ_reduction)%>%
  ggplot()+
  geom_bar(aes(x= community, y= mean, fill= community), color = "black", stat = "identity") +
  scale_fill_manual(labels= c("Mediterranean \n Alpine","Temperate \n Alpine"), values = c("darkgoldenrod1", "forestgreen"))+
  geom_errorbar(aes(x= community, y=mean, ymin=mean-se, ymax=mean+se), color = "black", linewidth=1, width = 0.4)+
  labs (y="Decreased germinationn (%)", subtitle= "b) Decreased germination by Habitat", x= "Habitat")+ #
  theme_classic (base_size = 9) + #theme_minimal for all species for mean treatment
  theme (axis.title.x = element_text(),
         axis.text.x = element_text (color="black"), #, angle = 20, vjust = 0.8
         legend.title = element_blank(),
         legend.position = "none")->fig3b;fig3b

### C) Water stress germination ordered by germination #####
finalgerm %>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  mutate(species = gsub("Helianthemum urrielense", "Helianthemum urrielense*", species))%>%
  filter(treatment== "C_alternate_WP"| treatment == "A_alternate_light")%>%
  select(community, species, habitat, community, treatment, opt_temp,petri,finalgerm, viable)%>%
  filter(!(species%in%nogerm_species$species))%>% # 7 species with 0 germ across all experiment
  filter(!(species%in%coldstrat_highgerm$species))%>% # 10 species with more than 50%germ in cold stratification
  filter (!species == "Solidago virgaurea")%>% # only germinated under cold stratification
  filter(!(species%in%nogerm_control_waterstress$species))%>%
  mutate(germpro = finalgerm/viable)%>% # generate Nan because of 0 viable
  mutate_all(~replace(., is.nan(.), 0))%>%
  select(community, species, habitat, community, treatment, petri,germpro)%>%
  spread(treatment, germpro)%>%
  mutate(germ_reduction=(1-(C_alternate_WP/A_alternate_light))*100)%>%
  mutate_all(~replace(., is.nan(.), 0))%>% # generate Nan because 0 germ pro
  convert_as_factor(species, community, petri, community)%>%
  mutate(community= fct_relevel(community, "Mediterranean", "Temperate" ))%>%
  mutate(community= fct_recode(community,"Mediterranean \n Alpine" ="Mediterranean", 
                               "Temperate \n Alpine"="Temperate" ))%>%
  mutate(species = ifelse(species=="Thymus praecox" & community == "Temperate \n Alpine", "Thymus praecox ", as.character(species)))%>%
  group_by(species, community)%>%
  get_summary_stats (germ_reduction)%>%
  ggplot(aes(x= reorder(species, mean), y= mean,fill = community)) +
  geom_bar (stat = "identity", color = "black") +
  geom_vline(xintercept = 40.5, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(labels= c("Mediterranean \n Alpine","Temperate \n Alpine"), values = c("darkgoldenrod1", "forestgreen"))+
  geom_errorbar( aes(species, mean, ymin = mean-se, ymax = mean+se), width = 0.5, linewidth = 1) +
  coord_flip() + 
  labs (subtitle = "c) Decreased germination by Species", y = "Decreased germination (%)")+
  theme_classic (base_size = 9) +
  theme (legend.position = "none",
         axis.title.y= element_blank(),
         axis.text.x= element_text(),
         axis.text.y = element_text(face= "italic"))-> fig3c;fig3c



# combine panels ###
library(patchwork)
((fig3a / fig3b) | fig3c)+
  plot_layout()-> fig3;fig3

ggsave(filename = "Fig4.png", plot =fig3 , path = "results/figures", 
       device = "png", dpi = 600, width = 155, height = 140, units = "mm") #, width = 180, units = "mm"

#ggpubr::ggarrange(fig2a, fig2b, ncol =2, nrow= 1,common.legend = FALSE, widths = c(1.5,1),align = "h")

