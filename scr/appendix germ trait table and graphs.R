library(ggplot2);library(tidyverse);library(rstatix);library(lubridate)
library(binom);library (vegan)

# Supporting information
# final germination percentage per species in each treatment
finalgerm %>%
  merge(read.csv("data/species.csv", sep=","), by= c("species", "code"))%>%
  group_by(species, community, treatment)%>%
  summarise(finalgerm=sum(finalgerm, na.rm=T),
            viable= sum(viable, na.rm=T))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  mutate(mean = round(mean, 2), 
         lower = round(lower, 2),
          upper = round(upper, 2))%>%
  mutate(speciesID = paste(species, community))%>%
  pivot_wider(id_cols = speciesID, names_from = treatment, values_from = c(n, mean, lower, upper), 
              names_glue = "{treatment}, {.value}")%>%
  separate_wider_delim(speciesID," ", names= c("genus", "species", "community"))%>%
  mutate(species= paste(genus, species))%>%
  dplyr::select(!genus)%>%
  dplyr::select(species, community, 
                `A_alternate_light, n`, `A_alternate_light, mean`,`A_alternate_light, lower`,`A_alternate_light, upper`,
                `B_alternate_dark, n`, `B_alternate_dark, mean`, `B_alternate_dark, lower`, `B_alternate_dark, upper`,
                `C_alternate_WP, n`,`C_alternate_WP, mean`,`C_alternate_WP, lower`,`C_alternate_WP, upper`,
                `D_constant_light, n`,`D_constant_light, mean`,`D_constant_light, lower`,`D_constant_light, upper`,
                `E_cold_stratification, n`,`E_cold_stratification, mean`,`E_cold_stratification, lower`,`E_cold_stratification, upper`)%>%
  write.csv("results/supplementary/Table S2. Species germination x treatment.csv", row.names = F)
  
# final germination percentage calculation per species and treatment ####
finalgerm%>%
  group_by(species, code, treatment)%>%
  summarise(finalgerm= sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code, treatment, finalgerm, viable,  mean, upper, lower)%>%
  as.data.frame()%>%
  mutate(treatment= as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light" , "Darkness"="B_alternate_dark", 
                               "Water stress"="C_alternate_WP", "Constant Temperatures"="D_constant_light", 
                               "Cold stratification" = "E_cold_stratification"))-> finalgerm_sp_treatment # final germ per petri dish 

# loop for individual graphs
for (var in unique(finalgerm_sp_treatment$species)) {
  treat_plot = ggplot(finalgerm_sp_treatment[finalgerm_sp_treatment$species==var,]) +
    geom_bar(aes(x = treatment, y = mean, ymin = 0, ymax = 1, fill = treatment), stat = "identity", color = "black") +
    geom_errorbar(aes(treatment, mean, ymin = lower, ymax = upper), width = 0.3, size =1.2) +
    scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature","Cold stratification" ), 
                      values = c("#2B83BA", "#ABDDA4", "#FFFFBF","#FDAE61", "#D7191C"))+
    labs(title= var, x = "Treatment", y = "Germination proportion") +
    geom_text(aes(x= treatment, y= -0.02, label=paste("n=", viable)),  size=4.5)+
    scale_y_continuous(limits = c(-0.02, 1),  breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    #facet_wrap(~ accession, nrow = 2) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text (size = 30, face= "italic"),
          #strip.text = element_text (size = 24, face = "italic"),
          axis.title.y = element_text (size=24), 
          axis.title.x = element_text (size=24), 
          axis.text.x= element_text (size=18, vjust = 0.5),
          legend.position = "none",
          plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm"))
  ggsave(treat_plot, file = paste0("Germination ", var,".png"),
         path = "results/supplementary/", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}

# separate individual plots for Thymus praecox
finalgerm_sp_treatment %>%
  filter(species=="Thymus praecox")%>%
  merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
  dplyr::select(community, species, treatment, mean, upper, lower, viable)%>%
  filter(community == "Temperate")%>%
  ggplot() +
  geom_bar(aes(x = treatment, y = mean, fill = treatment), stat = "identity", color = "black") +
  geom_errorbar(aes(treatment, mean, ymin = lower, ymax = upper), width = 0.3, size =1.2) +
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature", "Cold stratification"), 
                    values = c("#2B83BA", "#ABDDA4", "#FFFFBF","#FDAE61", "#D7191C"))+
  geom_text(aes(x= treatment, y= -0.02, label=paste("n=", viable)),  size=4.5)+ #position=position_stack(vjust=0.9),
  labs(title= "Temperate Thymus praecox", x = "Treatment", y = "Germination proportion") +
  scale_y_continuous(limits = c(-0.02, 1),  breaks = c(0, 0.25, 0.50, 0.75, 1)) +
  #facet_wrap(~ accession, nrow = 2) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 30),
        #strip.text = element_text (size = 24, face = "italic"),
        axis.title.y = element_text (size=24), 
        axis.title.x = element_text (size=24), 
        axis.text.x= element_text (size=18, vjust = 0.5),
        legend.position = "none",
        plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm"))->treat_plot;treat_plot
x11()
ggsave(treat_plot, file = "Germination Thymus praecox Temperate.png",
       path = "results/Supplementary", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
