library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)

# 1-check species names to match####
# dataframe with species data from germination experiments
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers)-> species 

# dataframe with species data from spatial survay
read.csv("data/spatial-survey-species-Med.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Med


read.csv("data/spatial-survey-species-Tem.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Tem

setdiff(spatial_sp_Med$species, species$species)
setdiff(spatial_sp_Tem$species, species$species)

# 2A- Cover/plot of species with germination traits ####
#MEDITERRANEAN
species %>%
  filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Med, by = "species")%>%
  filter (!species == "Solidago virgaurea")%>% # species with 0 germ
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover))%>%
  filter(rel_cover>79)%>%
  tally() 
# 15 plots with more than 80% coverage with sp traits (drivers) (considering raw cover, not to 100%)
# 65 plots with more than 80% coverage with sp traits (drivers) (considering relative cover, up to 100%)
# 64 plots with more than 80% coverage with sp traits (drivers + phenology)

#TEMPERATE
species %>%
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Tem, by = "species")%>%
  #filter (!species == "Euphrasia salisburgensis")%>%
  #filter (!species == "Gentiana verna")%>%
  #filter (!species == "Gentianella campestris")%>%
  #filter (!species == "Kobresia myosuroides")%>%
  #filter (!species == "Salix breviserrata")%>%
  #filter (!species == "Sedum album")%>%
  #filter (!species == "Sedum atratum")%>%
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover))%>%
  filter(rel_cover>79)%>%
  tally()
# 1 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)
# 48 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)
# 40 plot with more than 80% coverage with sp traits (drivers + phenology) (considering raw cover, not to 100%)

# 2B- N of species with germination traits per plot ####
#MEDITERRANEAN
species %>%
  filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Med, by = "species")%>%
  filter (!species == "Solidago virgaurea")%>% # species with 0 germ
  select(species, plot, N_sp)%>%
  group_by(plot) %>%
  summarise(N_sp = first(N_sp), 
            sp_traits = length(unique(species)), 
            sp_covered = sp_traits/N_sp)%>%
  filter(sp_covered>0.75) %>%
  print(n=84)
# in 65 plots more than 75% sp covered with drivers traits (focused on presence)
# in 57 plots more than 75% sp covered with drivers+phenology traits

#TEMPERATE
species %>%
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Tem, by = "species")%>%
  #filter (!species == "Euphrasia salisburgensis")%>%
  #filter (!species == "Gentiana verna")%>%
  #filter (!species == "Gentianella campestris")%>%
  #filter (!species == "Kobresia myosuroides")%>%
  #filter (!species == "Salix breviserrata")%>%
  #filter (!species == "Sedum album")%>%
  #filter (!species == "Sedum atratum")%>%
  select(species, plot, N_sp)%>%
  group_by(plot) %>%
  summarise(N_sp = first(N_sp), 
            sp_traits = length(unique(species)), 
            sp_covered = sp_traits/N_sp)%>%
  filter(sp_covered>0.75) %>%
  print(n=84)
# only 9 plots more than 75% sp covered with drivers traits (focused on presence)
# in 5 (2) plots more than 75% sp covered with drivers+phenology traits