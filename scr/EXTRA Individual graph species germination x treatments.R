library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix)

###UPDATE WITH ERROR BARS ####
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
  group_by(species, code, treatment)%>%
  summarise(finalgerm = sum(finalgerm),
            viable= sum(viable))%>%
  filter(!treatment=="E_cold_stratification")%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code, treatment, finalgerm, viable,  mean, upper, lower)%>%
  as.data.frame()-> finalgerm_sp_treatment # final germ per petri dish 


for (var in unique(finalgerm_sp_treatment$species)) {
  treat_plot = ggplot(finalgerm_sp_treatment[finalgerm_sp_treatment$species==var,], 
                      aes(x = treatment, y = mean, ymin = 0, ymax = 1, fill = treatment)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(treatment, mean, ymin = lower, ymax = upper), width = 0.3, size =1.2) +
    scale_fill_viridis_d() +
    labs(title= var, x = "Treatment", y = "Germination proportion") +
    scale_y_continuous(limits = c(0, 1),  breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    #facet_wrap(~ accession, nrow = 2) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text (size = 30),
          #strip.text = element_text (size = 24, face = "italic"),
          axis.title.y = element_text (size=24), 
          axis.title.x = element_text (size=24), 
          axis.text.x= element_text (size=18, vjust = 0.5),
          legend.position = "none",
          plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm"))
  ggsave(treat_plot, file = paste0("Germination ", var,".png"),
         path = NULL, scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}


# germination curves / rate  x species and x treatment #### 
#loop for many graphs   
viables_petri%>%
  group_by (species,  treatment) %>%
  summarise(viable = sum(viable))-> viable_sp

raw_df %>%
  mutate (D0 = 0) %>%
  select (species, treatment, initial, D0, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D0:D42) %>%
  #filter (!(treatment == "alternate_dark")) %>%
  mutate (scores = factor(scores, levels = c("D0","D7", "D14", "D21", "D28","D35", "D42"))) %>%
  arrange(species,  treatment, scores)%>% 
  group_by (species, treatment, scores)%>%
  mutate (germ = replace_na (germ, 0)) %>%
  summarise (germinated = sum(germ))%>%
  mutate (cum_germ = cumsum(germinated))  %>% #
  merge(viable_sp) %>%
  mutate(germpro = (cum_germ/viable), # 
         germpro = round (germpro, digit =2))%>% 
  as.data.frame()-> cumgerm_sp_treatment
cumgerm_sp_treatment$species <- as.character(cumgerm_sp_treatment$specie)
str(cumgerm_sp_treatment)

for (var in unique(cumgerm_sp_treatment$species)) {
  treat_plot = ggplot(cumgerm_sp_treatment[cumgerm_sp_treatment$species==var,], 
                      aes(x= scores, y= germpro, color= treatment, group = treatment))+
    geom_line(linewidth = 2) +
    ylim (0,1)+
    scale_color_manual (name= "Treatment", values = c ( "A_alternate_light" = "#440154FF", "B_alternate_dark" = "#2A788EFF","C_alternate_WP"= "#5ec962" ,"D_constant_light" = "#fde725")) +
    labs (title= var, x= "Germination scores", y ="Germination proportion")+
    theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
    theme (plot.title = element_text ( size = 32), #hjust = 0.5,
           legend.position = "right")
  theme(plot.title = element_text (size = 30),
        #strip.text = element_text (size = 24, face = "italic"),
        axis.title.y = element_text (size=24), 
        axis.title.x = element_text (size=24), 
        axis.text.x= element_text (size=18, vjust = 0.5),
        legend.position = "none",
        plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm"))
  ggsave(treat_plot, file = paste0("Germination ", var,".png"),
         path = NULL, scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}
