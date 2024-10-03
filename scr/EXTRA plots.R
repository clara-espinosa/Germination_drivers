library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix)

######## EXTRA PLOTS playing with VISUALIZATION #####
###UPDATE WITH ERROR BARS 
# germination curves / rate  x community and x treatment #### 
x11()
raw_df %>%
  mutate ( D0 = 0) %>%
  dplyr::select (species, treatment, initial, D0, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D0:D42) %>%
  filter (!(treatment == "alternate_dark")) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  mutate (scores = factor(scores, levels = c("D0","D7", "D14", "D21", "D28","D35", "D42"))) %>%
  arrange(species, treatment, scores)%>% 
  merge (species) %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianela campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")%>%
  filter (!species == "Teesdalia conferta")%>%
  group_by (community, treatment, scores)%>%
  summarise (germinated = sum(germ))%>%
  mutate (cum_germ = cumsum(germinated))  %>% #
  merge(viables_community) %>%
  mutate(germpro = (cum_germ/viable), # 
         germpro = round (germpro, digit =2)) %>%
  filter(!treatment == "B_alternate_dark")%>%
  ggplot(aes(x= scores, y= germpro, color= treatment, group = treatment))+
  geom_line(linewidth = 1.5) +
  facet_grid(~community) +
  #ylim (0,1)+
  scale_color_manual (name= "Treatment", values = c ("C_alternate_WP"= "#fde725", "A_alternate_light" ="#21918c", "D_constant_light" = "#5ec962")) +
  #scale_color_viridis_d() +
  labs (title = "Cumulative Germination")+
  xlab ("Germination scores") +
  ylab("Germination proportion") +
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         legend.position = "right")

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
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code, treatment, finalgerm, viable,  mean, upper, lower)%>%
  as.data.frame()%>%
  mutate(treatment= as.factor(treatment))%>%
  mutate(treatment= fct_recode(treatment, "Control"="A_alternate_light" , "Darkness"="B_alternate_dark", 
                               "Water stress"="C_alternate_WP", "Constant Temperatures"="D_constant_light"))-> finalgerm_sp_treatment # final germ per petri dish 


for (var in unique(finalgerm_sp_treatment$species)) {
  treat_plot = ggplot(finalgerm_sp_treatment[finalgerm_sp_treatment$species==var,]) +
    geom_bar(aes(x = treatment, y = mean, ymin = 0, ymax = 1, fill = treatment), stat = "identity", color = "black") +
    geom_errorbar(aes(treatment, mean, ymin = lower, ymax = upper), width = 0.3, size =1.2) +
    scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                      values = c("#2A788EFF", "#440154FF", "#FDE725FF","#35B779FF"))+
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


# germination curves / rate  x species and x treatment IF USED NEEDS TO BE UPDATED #### 
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

# cold stratification data visualization ####
x11()
cold_strat %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  dplyr::select (species, code,  treatment, finalgerm, viable, mean, upper, lower)%>%
  group_by (species)%>%
  summarize(germpro = mean(mean), 
            upper = mean(upper),
            lower = mean (lower)) %>% 
  ggplot(aes(x= species, y= germpro, fill = species)) +
  geom_bar (stat = "identity") +
  geom_errorbar( aes(species, germpro, ymin = lower, ymax = upper), width = 0.5, linewidth = 0.5) +
  coord_flip() + 
  labs (title = "Germination during cold stratification", x = "Species", y = "Germination proportion")+
  theme_minimal (base_size = 16) +
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         legend.position = "none")


# GDD/FDD correlation with treatment germination responses #####
##GDD        
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>% #-> mcmc_ibc
  merge(species2)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  ggplot(aes(x= GDD, y= mean, fill= treatment), color = "black")+
  geom_point(stat = "identity", shape=21, size=4) +
  facet_grid(~treatment) + 
  geom_smooth(method = "loess")+
  ylim (-0.1,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "dimgrey", "chocolate1","#7AD151FF")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "GDD (Growing dregree days)", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =16),
         axis.title.x= element_text (size =18), 
         axis.text.x= element_blank(), 
         legend.position = "bottom")
##FDD
finalgerm %>%
  group_by (species, code, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  arrange(species, treatment) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viable_sp, by= c("code", "species", "treatment"))%>% #-> mcmc_ibc
  merge(species2)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  ggplot(aes(x= FDD, y= mean, fill= treatment), color = "black")+
  geom_point(stat = "identity", shape=21, size=4) +
  facet_grid(~treatment) +
  geom_smooth(method = "lm")+
  ylim (-0.1,0.5)+
  scale_fill_manual(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature"), 
                    values = c("#2A788EFF", "dimgrey", "chocolate1","#7AD151FF")) + 
  #"Darkness" = "dimgrey", "Water stress"= "chocolate1", "Constant Temperature" = "#7AD151FF", "Cold stratification"
  labs (x = "FDD (Freezing degree days)", y="Germination proportion")+
  theme_classic (base_size = 18) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =16),
         axis.title.x= element_text (size =18), 
         axis.text.x= element_blank(), 
         legend.position = "none")

