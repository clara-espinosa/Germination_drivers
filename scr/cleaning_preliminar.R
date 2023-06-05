library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis)

# dataframe with raw data ####
raw_df <- read_excel("C:/Users/Usuario/OneDrive - Universidad de Oviedo/IMIB/Experiments/Light_Dark_Temperature/raw_data.xlsx") 
raw_df$species <- factor(raw_df$species) 
raw_df$code <- factor(raw_df$code) 
raw_df$treatment <- factor(raw_df$treatment) 
raw_df$temp <- factor(raw_df$temp) 
str(raw_df)

# dataframe with species data ####
species <- read_excel("C:/Users/Usuario/OneDrive - Universidad de Oviedo/IMIB/Experiments/Light_Dark_Temperature/species.xlsx")
species %>%
  select(species, temperature_regime, family, community, habitat) -> species
str(species)

# species with 0 germination across all experiment
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianela campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix fontqueri")%>%
  filter (!species == "Sedum album CF")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")
  
# viables calculation x petri dish ####
raw_df %>%
  select (species, treatment, petri, initial, empty, fungus) %>%
  group_by (species, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  mutate (viablePER = (viable/initial)*100) %>%
  select (species, treatment, petri, viable, viablePER) -> viables_petri

# D0 germination check change to cold_stratification treatment ####
raw_df %>%
  merge(viables_petri) %>%
  select (species, petri, viable, initial, cold_strat) %>% ## empty, fungus to substract to initial
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species) %>%
  summarize(germ_under_snow = mean(cold_strat/viable)) %>% 
  mutate (germ_under_snow= replace_na(germ_under_snow, 0)) %>% # salix and solidago have some petridish with all fungus
  arrange (germ_under_snow) %>%
  mutate (germpro = germ_under_snow, 
        treatment = "cold_stratification") %>%
  mutate(treatment = as.factor(treatment)) %>%
  select (species,  treatment, germpro)-> cold_strat

x11()
cold_strat %>%
  group_by (species)%>%
  summarize(germpro = mean(germpro)) %>% 
  ggplot(aes(x= species, y= germpro, fill = species)) +
  geom_bar (stat = "identity") +
  coord_flip() + 
  labs (title = "Germination during cold stratification")+
  xlab ("Species") +
  ylab("Germination proportion") +
  theme_minimal (base_size = 16) +
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
    legend.position = "none")

# final germination percentage calculation discounting cold stratification germ ####
raw_df %>%
  select (species, treatment, petri, initial, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species,  treatment, petri)%>%
  summarize (final_cum_germ = sum(germ), 
             initial = first(initial)) %>% # 
  merge(viables_petri) %>%
  mutate(germpro = (final_cum_germ/viable), # 
         germpro = round (germpro, digit =2)) %>%
  select (species, treatment, petri, final_cum_germ, viable, germpro)-> finalgerm # final germ per petri dish 

# graph comparing treatments !!!! ####
finalgerm %>%
  group_by (species, treatment) %>%
  summarize(germpro = mean (germpro)) %>% # mean to match dataset and substract cold stratification germ
  rbind(cold_strat) %>%
  spread(treatment, germpro)%>%
  mutate( alternate_dark = alternate_dark - cold_stratification)  %>% # substract cold stratification
  mutate (alternate_dark = ifelse(alternate_dark < 0, 0, alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "germpro", alternate_dark:cold_stratification) %>%
  merge (species, by = "species") %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianela campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix fontqueri")%>%
  filter (!species == "Sedum album CF")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  na.omit() %>%
  group_by (treatment, community) %>% #, temperature_regime 
  summarize (germpro = mean(germpro)) %>%
  mutate (treatment = factor(treatment, 
                             levels = c("cold_stratification", "alternate_dark", "alternate_light", 
                                        "constant_light", "alternate_WP" ))) %>%
  ggplot()+
  geom_bar(aes(x= treatment, y= germpro, fill= treatment), stat = "identity") +
  facet_grid(~community) +
  ylim (0,1)+
  scale_fill_viridis_d() +
  labs (title = "Mean Germination per Treatment")+
  xlab ("Treatment") +
  ylab("Germination proportion") +
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         legend.position = "none")

#loop for many graphs 
finalgerm %>%
  group_by (species, treatment) %>%
  summarize(germpro = mean (germpro)) %>% # mean to match dataset and substract cold stratification germ
  rbind(cold_strat) %>%
  spread(treatment, germpro)%>%
  mutate( alternate_dark = alternate_dark - cold_stratification)  %>%
  mutate (alternate_dark = ifelse(alternate_dark < 0, 0, alternate_dark)) %>%
  gather ("treatment", "germpro", alternate_dark:cold_stratification)%>%
  mutate (germpro = replace_na (germpro, 0))%>%
  group_by (species, treatment) %>%
  summarize ( germpro = mean(germpro))%>%
  as.data.frame()-> finalgerm_sp_treatment
str(finalgerm_sp_treatment)
finalgerm_sp_treatment$species <- as.character(finalgerm_sp_treatment$specie)

for (var in unique(finalgerm_sp_treatment$species)) {
  treat_plot = ggplot(finalgerm_sp_treatment[finalgerm_sp_treatment$species==var,], 
                      aes(x = treatment, y = germpro, ymin = 0, ymax = 1, fill = treatment)) +
    geom_bar(stat = "identity") +
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
viables_petri %>%
  merge(species)%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianela campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix fontqueri")%>%
  filter (!species == "Sedum album CF")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")%>%
  group_by (community, treatment)%>%
  summarise(viable= sum(viable))-> viables_community
x11()
raw_df %>%
  mutate ( D0 = 0) %>%
  select (species, treatment, initial, D0, D7, D14, D21, D28, D35, D42) %>%
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
  filter (!species == "Salix fontqueri")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")%>%
  group_by (community, treatment, scores)%>%
  summarise (germinated = sum(germ))%>%
  mutate (cum_germ = cumsum(germinated))  %>% #
  merge(viables_community) %>%
  mutate(germpro = (cum_germ/viable), # 
         germpro = round (germpro, digit =2)) %>%
  ggplot(aes(x= scores, y= germpro, color= treatment, group = treatment))+
  geom_line(size = 1.5) +
  facet_grid(~community) +
  ylim (0,1)+
  scale_color_manual (name= "Treatment", values = c ("alternate_WP"= "#fde725", "alternate_light" ="#21918c", "constant_light" = "#5ec962")) +
  #scale_color_viridis_d() +
  labs (title = "Cumulative Germination")+
  xlab ("Germination scores") +
  ylab("Germination proportion") +
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         legend.position = "right")
x11()
#loop for many graphs   
viables_petri%>%
  group_by (species,  treatment) %>%
  summarise(viable = sum(viable))-> viable_sp

raw_df %>%
  mutate (D0 = 0) %>%
  select (species, treatment, initial, D0, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D0:D42) %>%
  filter (!(treatment == "alternate_dark")) %>%
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
    geom_line(size = 1.5) +
    ylim (0,1)+
    scale_color_manual (name= "Treatment", values = c ("alternate_WP"= "#fde725", "alternate_light" ="#21918c", "constant_light" = "#5ec962")) +
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
