library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(binom)

# dataframe with raw data ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) -> raw_df
str(raw_df)
unique(raw_df$species)

# dataframe with species data ####
read.csv("data/species.csv", sep=";") %>%
  select(species, code, temperature_regime, family, community, habitat)  %>%
  convert_as_factor(species, code, temperature_regime, family, community, habitat)  -> species
str(species)
unique(species$family)

# species with 0 germination across all experiment (remove from analysis)
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")
  
# viables calculation x petri dish ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, empty, fungus) %>%
  group_by (species, code, treatment, petri) %>%
  mutate (viable = initial -(empty + fungus)) %>%
  mutate (viablePER = (viable/initial)*100) %>%
  select (species, code, treatment, petri, viable, viablePER) -> viables_petri
# viables x community x treatment
  read.csv("data/raw_data.csv", sep = ";") %>%
    convert_as_factor(species, code, treatment, temp) %>%
    select (species, code, treatment, petri, initial, empty, fungus) %>%
    group_by (species, code, treatment, petri) %>%
    mutate (viable = initial -(empty + fungus)) %>%
    mutate (viablePER = (viable/initial)*100) %>%
    select (species, code, treatment, petri, viable, viablePER)%>%
    merge(species)%>%
    group_by (community, treatment)%>%
    summarise(viable= sum(viable))-> viables_community
# viables x species x treatment
read.csv("data/raw_data.csv", sep = ";") %>%
    convert_as_factor(species, code, treatment, temp) %>%
    select (species, code, treatment, petri, initial, empty, fungus) %>%
    group_by (species, code, treatment, petri) %>%
    mutate (viable = initial -(empty + fungus)) %>%
    mutate (viablePER = (viable/initial)*100) %>%
    select (species, code, treatment, petri, viable)%>%
    rbind (cold_strat_viables)%>%
    merge(species)%>%
    group_by (community, species,treatment)%>%
    summarise(viable= sum(viable))-> viable_sp

# list of accession with less than 25% of viable seeds (remove from analysis) 
viables_petri %>%
  group_by(species, code, treatment) %>%
  summarise(viable = sum(viable), viablePER = mean(viablePER))%>%
  filter(viablePER<25) # salix breviserrata


# D0 germination check change to cold_stratification treatment and mean values ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  merge(viables_petri) %>%
  select (species, code, petri, viable, initial, cold_strat) %>% ## 
  na.omit() %>% #remove alternate dark treatment rows
  group_by (species, code,petri) %>%
  summarize(finalgerm = sum(cold_strat),
            viable = sum(viable)) %>% 
  #mutate (viablePER = (viable/initial)*100) %>%
  mutate (finalgerm = replace_na(finalgerm , 0)) %>% # salix and solidago have some petridish with all fungus
  mutate (treatment = "E_cold_stratification") %>%
  #select(species, code, treatment, petri, viable)->cold_strat_viables
  mutate(treatment = as.factor(treatment)) %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code,  treatment, finalgerm, viable, mean, upper, lower)-> cold_strat

x11()
cold_strat %>%
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

# WORK on IT !!! final germination percentage calculation discounting cold stratification germ ####
read.csv("data/raw_data.csv", sep = ";") %>%
  convert_as_factor(species, code, treatment, temp) %>%
  select (species, code, treatment, petri, initial, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D7:D42) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  group_by (species, code, treatment, petri)%>%
  summarize (finalgerm = sum(germ), 
             initial = first(initial)) %>% # 
  merge(viables_petri) %>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  select (species, code, treatment, petri, finalgerm, viable, viablePER, mean, upper, lower)-> finalgerm # final germ per petri dish 

# graph comparing final germination between treatments !!!! ####
### substract cold_strat from dark treatment, and add binomial errors!!
cold_strat%>%
  merge (species)%>%
  group_by (community, species, treatment) %>%
  summarize(finalgerm = sum (finalgerm))-> cold_strat2

finalgerm %>%
  merge(species)%>%
  group_by (community, species, treatment) %>%
  summarize(finalgerm = sum (finalgerm)) %>% # mean to match dataset and substract cold stratification germ
  rbind(cold_strat2)%>%
  arrange(species, treatment) %>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  merge(viable_sp, by= c("community", "species", "treatment"))%>%
  #merge (species, by = "species") %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  na.omit() %>%
  #group_by (treatment, community) %>% #, temperature_regime 
  summarise(finalgerm=sum(finalgerm),
            viable = sum(viable))%>%
  mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
  ggplot()+
  geom_bar(aes(x= treatment, y= mean, fill= treatment), stat = "identity") +
  geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", size=1, width = 0.4)+
  facet_grid(~community) +
  #ylim (0,1)+
  scale_fill_viridis_d(name = "Treatments",labels = c("Control", "Darkness", "Water stress", "Constant Temperature", "Cold stratification")) +
  labs (title = "Mean Germination per Treatment", x = "Treatments", y="Germination proportion")+
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         strip.text = element_text (size =20),
         axis.text.x= element_blank(), 
         legend.position = "right")

# germination curves / rate  x community and x treatment #### 
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
  filter(!treatment == "B_alternate_dark")%>%
  ggplot(aes(x= scores, y= germpro, color= treatment, group = treatment))+
  geom_line(linewidth = 1.5) +
  facet_grid(~community) +
  ylim (0,1)+
  scale_color_manual (name= "Treatment", values = c ("C_alternate_WP"= "#fde725", "A_alternate_light" ="#21918c", "D_constant_light" = "#5ec962")) +
  #scale_color_viridis_d() +
  labs (title = "Cumulative Germination")+
  xlab ("Germination scores") +
  ylab("Germination proportion") +
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         legend.position = "right")
