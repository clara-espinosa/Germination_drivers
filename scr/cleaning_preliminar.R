library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix)

# dataframe with raw data ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) -> raw_df
str(raw_df)

# dataframe with species data ####
read.csv("data/species.csv", sep=";") %>%
  select(species, temperature_regime, family, community, habitat)  %>%
  convert_as_factor(species, temperature_regime, family, community, habitat)  -> species
str(species)
unique(species$family)
# species with 0 germination across all experiment (remove from analysis)
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

  # list of accesion with less than 25% of viable seeds (remove from analysis) 
viables_petri %>%
  group_by(species, treatment) %>%
  summarise(viable = sum(viable), viablePER = mean(viablePER))%>%
  filter(viablePER>25)

viables_petri%>%
  group_by(species)%>%
  summarise(viable=sum(viable))%>%
  mutate(total_viable= sum(viable))

# D0 germination check change to cold_stratification treatment ####
raw_df %>%
  merge(viables_petri) %>%
  select (species, petri, viable, initial, cold_strat) %>% ## empty, fungus to subtract to initial
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

# graph comparing final germination between treatments !!!! ####
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

# germination curves / rate  x species and x treatment #### 
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
