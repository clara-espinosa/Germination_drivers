library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(lubridate);library(vegan)
library(rstatix);library(geomtextpath);library(psych)
library (ggrepel);library (ggpubr) ;library(binom)
library(RColorBrewer);library(ggtext);library(patchwork)

# cold stratification visualization
### A) cold stratification ordered by germination #####
cold_strat %>%
  merge (species)%>%
  filter(!(species%in%nogerm_species$species))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  dplyr::select (community, family, species, treatment, biogeography, finalgerm, viable, mean, upper, lower)%>%
  group_by (biogeography, family, species)%>%
  summarize(germpro = mean(mean), 
            upper = mean(upper),
            lower = mean (lower)) %>% 
  #mutate(cold= ifelse(germpro>0.1, "yes", "no"))%>% 
  #ggplot(aes(x= new_species, y= germpro,fill = cold)) +
  ggplot(aes(x= reorder(species, germpro), y= germpro,fill = biogeography)) +
  geom_bar (stat = "identity", color = "black") +
  geom_vline(xintercept = 40.5, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("Mediterranean"="#FDAE61","Broad range"= "#FFFFBF","Endemic"= "#ABDDA4", "Eurosiberian"= "#2B83BA"))+
  geom_errorbar( aes(species, germpro, ymin = lower, ymax = upper), width = 0.5, linewidth = 0.5) +
  annotate("segment", x = 41, y = 0.75, xend = 45, yend = 0.75, arrow = arrow(type = "closed", length = unit(0.35, "cm")))+
  annotate("label",x=43,y=0.89, label = "Under snow \n germination", size= 3)+
  annotate("segment", x = 40, y = 0.75, xend = 36, yend = 0.75, arrow = arrow(type = "closed", length = unit(0.35, "cm")))+
  annotate("label",x=38,y=0.89, label = "Spring \n germination", size= 3)+
  coord_flip() + 
  labs (title = "Germination during cold stratification", subtitle = "A) Individual species", y = "Germination proportion")+
  theme_classic (base_size = 10) +
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         legend.position = "none",
         axis.title.y= element_blank(),
         axis.text.x= element_text(),
         axis.text.y = element_text(face= "italic", size = 9))-> fig2a;fig2a
### B) differences between biogeography ####
x11()
cold_strat%>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm), viable = sum(viable)) %>%
  merge (species)%>%
  filter(!(species%in%nogerm_species$species))%>%
  group_by(biogeography)%>%
  summarize(finalgerm = sum (finalgerm), 
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  merge(mcmc_cold%>%
          group_by (species,treatment) %>%
          summarize(finalgerm = sum (finalgerm), viable = sum(viable)) %>%
          filter(!(species%in%nogerm_species$species))%>% 
          merge (species)%>%
          group_by(biogeography)%>%
  tally(), by = "biogeography")%>%
  rename(n.seeds = n.x,
         n.species= n.y)%>%
  convert_as_factor(biogeography)%>%
  mutate(biogeography= fct_relevel(biogeography, "Mediterranean", "Eurosiberian", "Endemic", "Broad range" ))%>%
  ggplot()+
  geom_bar(aes(x= biogeography, y= mean, fill= biogeography), color = "black", stat = "identity") +
  scale_fill_manual(values = c("Mediterranean"="#FDAE61","Broad range"= "#FFFFBF","Endemic"= "#ABDDA4", "Eurosiberian"= "#2B83BA"))+
  geom_errorbar(aes(x= biogeography, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4)+
  geom_segment (aes(x= 1,xend =2,  y = 0.3, yend= 0.3), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.31, label = "*", size= 6)+
  geom_segment (aes(x= 1,xend =3,  y = 0.33, yend= 0.33), color = "black", linewidth = 1.3, show.legend = F)+
  annotate ("text", x= 2, y= 0.34, label = "**", size= 6)+
  geom_text(aes(x= biogeography, y= 0.01, label=paste("n =", n.species)),  size=3)+
  ylim (0,0.35)+
  labs (y="Germination proportion", subtitle= "B) Biogeographical distribution")+ #
  theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         axis.title.x = element_blank(),
         axis.text.x = element_text (color="black", angle = 20, vjust = 0.8),
         legend.title = element_blank(),
         legend.position = "none")-> fig2b;fig2b

# combine panels ###

fig2a+ fig2b + 
  plot_layout(widths = c(1.5,1))-> fig2;fig2

ggsave(filename = "cold stratification.png", plot =fig2 , path = "results/figures", 
       device = "png", dpi = 600) #, width = 180, units = "mm"

#ggpubr::ggarrange(fig2a, fig2b, ncol =2, nrow= 1,common.legend = FALSE, widths = c(1.5,1),align = "h")

