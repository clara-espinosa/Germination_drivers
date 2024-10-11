library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(vegan);library(geomtextpath)

# visualization exploration with density plots
##### MEDITERRANEAN SPECIES ALONG MICROCLIMATIC GRADIENTS ###########
# dataframe with species data from germination experiments
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers) %>%
  filter (community == "Mediterranean")%>%
  #filter(!species=="Solidago virgaurea")%>% # 0 germ across all treatments
  #filter(!species == "Teesdalia conferta") %>% # all germinated during cold stratification
  #filter(!species == "Jasione laevis")%>% # without leaves traits
  as.data.frame()-> species_med

# dataframe with species data from spatial survay
read.csv("data/spatial-survey-species-Med.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Med

# 3-  plot x Environmental data
read.csv("data/spatial-survey-temperatures-Med.csv") %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M")) %>%
  group_by(plot, Day = lubridate::floor_date(Time, "day")) %>%
  dplyr::summarize(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 5, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  group_by(plot, Month = lubridate::floor_date(Day, "month")) %>%
  dplyr::summarize(T = mean(T), X = mean(X), N = mean(N), # Daily mean, max, min
                   Snow = sum(Snow), # Snow days per month
                   FreezeThaw = sum(FreezeThaw), # Freeze-thaw days per month
                   FDD = sum(FDD), # FDD per month
                   GDD = sum(GDD)) %>% # GDD per month
  group_by(plot) %>%
  dplyr::summarize(bio1 = round(mean(T),2), # Annual Mean Temperature
                   bio2 = round(mean(X - N),2), # Mean Diurnal Range (Mean of monthly (max temp - min temp))
                   bio7 = round(max(X) - min(N), 2), # Temperature Annual Range (BIO5-BIO6)
                   Snw = sum(Snow),
                   FDD = round(abs(sum(FDD)),2), # FDD per year
                   GDD = round(sum(GDD)),2) %>%  # GDD per year
  merge(read.csv("data/spatial-survey-header-Med.csv")) %>% 
  dplyr::select(plot, elevation,Snw:GDD)  %>% 
  merge(spatial_sp_Med, by = "plot") %>% 
  dplyr::select(plot, elevation,Snw:GDD, species, cover, rel_cover)%>% 
  merge(sp_med, by="species")%>%
  dplyr::select(plot, elevation,Snw:GDD, species, cover, rel_cover, family)%>%
  gather(env, micro_value, elevation:GDD)%>%
  ggplot()+
  geom_density(aes(x=micro_value, group=family, fill=family), alpha = 0.5)+ #, position = "fill"
  geom_textdensity (aes(x=micro_value, group=family, label = family), hjust = "ymax", position= "jitter", color="black", show.legend = F)+
  #geom_density(aes(x=micro_value, group=species, fill=species), alpha = 0.5)+ #, position = "fill"
  #geom_textdensity (aes(x=micro_value, group=species, label = species), hjust = "ymax", position= "jitter", color="black", show.legend = F)+
  facet_wrap(~env, scales = "free", ncol=1)+
  theme_classic()
x11()

##### MEDITERRANEAN PLOTS ALONG MICROCLIMATIC GRADIENTS ###########
# join traits subset
germ_odds_M%>% # log odds ratios
  merge(seed_mass_M)%>% # log transformed
  merge(plant_height_M)%>% # log transformed
  merge(leaves_traits_M )%>% #only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height,leaf_area, LDMC, SLA,
                odds_B_dark, odds_C_WP, odds_D_constant)-> species_traits_M
# 3-  plot x Environmental data 
read.csv("data/spatial-survey-temperatures-Med.csv") %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M")) %>%
  group_by(plot, Day = lubridate::floor_date(Time, "day")) %>%
  dplyr::summarize(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 5, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  group_by(plot, Month = lubridate::floor_date(Day, "month")) %>%
  dplyr::summarize(T = mean(T), X = mean(X), N = mean(N), # Daily mean, max, min
                   Snow = sum(Snow), # Snow days per month
                   FreezeThaw = sum(FreezeThaw), # Freeze-thaw days per month
                   FDD = sum(FDD), # FDD per month
                   GDD = sum(GDD)) %>% # GDD per month
  group_by(plot) %>%
  dplyr::summarize(bio1 = round(mean(T),2), # Annual Mean Temperature
                   bio2 = round(mean(X - N),2), # Mean Diurnal Range (Mean of monthly (max temp - min temp))
                   bio7 = round(max(X) - min(N), 2), # Temperature Annual Range (BIO5-BIO6)
                   Snw = sum(Snow),
                   FDD = round(abs(sum(FDD)),2), # FDD per year
                   GDD = round(sum(GDD)),2) %>%  # GDD per year
  merge(read.csv("data/spatial-survey-header-Med.csv")) %>% 
  #dplyr::select(plot, site, aspect, elevation, latitude, longitude, bio1:GDD)%>% 
  dplyr::select(site,plot, elevation,bio1:GDD)%>% 
  filter(plot%in%spatial_env_med$plot)%>%
  merge(spatial_sp_Med, by = "plot") %>% 
  dplyr::select(site,plot, elevation,Snw:GDD, species, cover, rel_cover)%>% 
  merge(sp_med, by="species")%>%
  merge(species_traits_M, by="species")%>%
  dplyr::select (site,plot, elevation:GDD, species, cover, rel_cover, seed_mass:odds_D_constant) %>%
  gather (traits, trait_value, seed_mass:odds_D_constant )%>%
  gather(env, micro_value, elevation:GDD)%>%
  ggplot()+
  geom_density(aes(x=trait_value, group=plot, fill=site), alpha = 0.5)+ #, position = "fill"
  geom_textdensity (aes(x=trait_value, group=plot, label = plot), hjust = "ymax", position= "jitter", color="black", show.legend = F)+
  #geom_density(aes(x=micro_value, group=species, fill=species), alpha = 0.5)+ #, position = "fill"
  #geom_textdensity (aes(x=micro_value, group=species, label = species), hjust = "ymax", position= "jitter", color="black", show.legend = F)+
  facet_wrap(~traits, scales = "free", ncol=1)+
  theme_classic()
