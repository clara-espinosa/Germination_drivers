library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(lubridate);library(vegan)
library(rstatix);library(geomtextpath);library(psych)
library (ggrepel);library(binom)
library(RColorBrewer);library(ggtext);library(patchwork)
library(scales);library(ggtext)

# cold stratification visualization
### A) cold stratification ordered by germination #####
cold_strat %>%
  merge (species)%>%
  mutate(species = gsub("Helianthemum urrielense", "Helianthemum urrielense*", species))%>%
  filter(!(species%in%nogerm_species$species))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  dplyr::select (community, family, species, treatment, community, finalgerm, viable, mean, upper, lower)%>%
  group_by (community, family, species)%>%
  summarize(germpro = mean(mean), 
            upper = mean(upper),
            lower = mean (lower)) %>% 
  #mutate(cold= ifelse(germpro>0.1, "yes", "no"))%>% 
  #ggplot(aes(x= new_species, y= germpro,fill = cold)) +
  ggplot(aes(x= reorder(species, germpro), y= germpro,fill = community)) +
  geom_bar (stat = "identity", color = "black") +
  geom_vline(xintercept = 40.5, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("Mediterranean"="darkgoldenrod1","Temperate"= "forestgreen"))+
  geom_errorbar( aes(species, germpro, ymin = lower, ymax = upper, color= community), width = 0.8, linewidth = 0.8) + #
  scale_color_manual(values = c("Mediterranean"="darkgoldenrod1","Temperate"= "forestgreen"))+
  annotate("segment", x = 41, y = 0.73, xend = 45, yend = 0.73, arrow = arrow(type = "closed", length = unit(0.28, "cm")))+
  annotate("label",x=43,y=0.9, label = "Under snow \n germination", size= 2.7)+
  annotate("segment", x = 40, y = 0.73, xend = 36, yend = 0.73, arrow = arrow(type = "closed", length = unit(0.28, "cm")))+
  annotate("label",x=38,y=0.9, label = "Spring \n germination", size= 2.7)+
  scale_x_discrete(labels = label)+
  coord_flip() + 
  labs (subtitle = "A) By Species", y = "Germination proportion")+
  theme_classic (base_size = 10) +
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         legend.position = "none",
         axis.title.y= element_blank(),
         axis.text.x= element_text(),
         axis.text.y = element_markdown(face= "italic", size = 9))-> fig2a;fig2a
### B) differences between community ####
x11()
cold_strat%>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm), viable = sum(viable)) %>%
  merge (species)%>%
  filter(!(species%in%nogerm_species$species))%>%
  group_by(community)%>%
  summarize(finalgerm = sum (finalgerm), 
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  merge(mcmc_cold%>%
          group_by (community, species) %>%
          summarize(finalgerm = sum (finalgerm), viable = sum(viable)) %>%
          group_by(community)%>%
  tally(), by = "community")%>%
  rename(n.seeds = n.x,
         n.species= n.y)%>%
  convert_as_factor(community)%>%
  mutate(community= fct_relevel(community, "Mediterranean", "Temperate" ))%>%
  mutate(community= fct_recode(community,"Mediterranean \n Alpine" ="Mediterranean", 
                               "Temperate \n Alpine"="Temperate" ))%>%
  ggplot()+
  geom_bar(aes(x= community, y= mean, fill= community), color = "black", stat = "identity") +
  scale_fill_manual(labels= c("Mediterranean \n Alpine","Temperate \n Alpine"), values = c("darkgoldenrod1", "forestgreen"))+
  geom_errorbar(aes(x= community, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+
  geom_segment (aes(x= 1,xend =2,  y = 0.27, yend= 0.27), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.28, label = "**", size= 6)+
  geom_text(aes(x= community, y= 0.01, label=paste("n =", n.species)),  size=3)+
  ylim (0,0.3)+
  labs (y="Germination proportion", subtitle= "B) By Habitat")+ #
  theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         axis.title.x = element_blank(),
         axis.text.x = element_text (color="black"), #, angle = 20, vjust = 0.7
         legend.title = element_blank(),
         legend.position = "none")-> fig2b;fig2b

# combine panels ###
library(patchwork)
fig2a+ fig2b + 
  plot_layout(widths = c(1.5,1))+ plot_annotation (title = "Germination during cold stratification (n=50)")-> fig2;fig2

ggsave(filename = "cold stratification.png", plot =fig2 , path = "results/figures", 
       device = "png", dpi = 600, width = 170, height = 150, units = "mm") #

#ggpubr::ggarrange(fig2a, fig2b, ncol =2, nrow= 1,common.legend = FALSE, widths = c(1.5,1),align = "h")

