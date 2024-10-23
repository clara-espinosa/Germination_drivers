library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(lubridate);library(vegan)
library(rstatix);library(geomtextpath);library(psych)
library (ggrepel);library (ggpubr)

######## EXTRA PLOTS playing with VISUALIZATION #####
# germination curves / rate  x community and x treatment (update with error bars) #### 
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

#final germination percentage calculation (update with error bars) ####
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


# visualization exploration with density plots #####
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

# loop for visualization plots CM x microclimatic gradients####
CM_M%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_M2)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CM_microclima_M

for (var in unique(CM_microclima_M$micro_variable)){
  micro_plot=ggplot(CM_microclima_M[CM_microclima_M$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Means in Mediterranean system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CM Mediterranean x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}
## facet for trait
x11()
for (var in unique(CM_microclima_M$trait)){
  micro_plot=ggplot(CM_microclima_M[CM_microclima_M$trait==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~micro_variable, ncol = 4, scales = "free")+
    labs(title="Community Means in Mediterranean system", subtitle = var, y= var )+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          axis.title.x = element_blank(),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CM Mediterranean x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 125, units = "mm", dpi = 600)
}

for (var in unique(CWM_microclima_M$micro_variable)){
  micro_plot=ggplot(CWM_microclima_M[CWM_microclima_M$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Weigthed Means in Mediterranean system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CWM Mediterranean x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}
for (var in unique(CM_microclima_T$micro_variable)){
  micro_plot=ggplot(CM_microclima_T[CM_microclima_T$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Means in Temperate system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CM Temperate x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}
for (var in unique(CWM_microclima_T$micro_variable)){
  micro_plot=ggplot(CWM_microclima_T[CWM_microclima_T$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Weighted Means in Temperate system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CWM Temperate x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}

# Visualization as de Bello 2013 foe community metrics (problem with FD, different scale) #####
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                 "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")

CM_microclima_M%>%
  rename(CM = value)%>%
  merge(CWM_microclima_M, by = c("site", "plot", "micro_variable", "value_micro", "trait"))%>%
  rename(CWM = value)%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  ggplot()+
  #geom_point(aes(y=CM, x= value_micro, color= "black" ),size =1, show.legend = T, alpha = 0.5)+
  geom_smooth (aes(y=CM, x= value_micro, color = "black"),method = "lm", se=F, inherit.aes=F, linewidth =1.5)+
  #geom_point(aes(y=CWM, x= value_micro,color= "grey"),  size =1, show.legend = T, alpha = 0.5)+
  geom_smooth (aes(y=CWM, x= value_micro, color = "grey"),method = "lm", se=F, inherit.aes=F, linewidth =1.5)+
  facet_grid(factor(trait, levels = c("odds_B_dark","odds_C_WP","odds_D_constant",
                                      "seed_mass","plant_height", 
                                      "leaf_area", "LDMC", "SLA"))~micro_variable, scales = "free", 
             labeller = as_labeller(trait_names))+
  geom_text(data=ann.sig.CM.M, label =ann.sig.CM.M$label, aes(x=x, y=y),size= 6, color = "black")+
  geom_text(data=ann.sig.CWM.M, label =ann.sig.CWM.M$label, aes(x=x, y=y),size= 6, color = "grey")+
  #labs(title="Community Means in Mediterranean system")+
  scale_color_manual (name = "Community metrics", values = c("black"="black", "grey"="grey"), 
                      labels = c("CM", "CWM"), guide= guide_legend(theme=theme(legend.title.position = "top")))+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        strip.background = element_blank(),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_blank(),
        #axis.text= element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(plot_CM_M, file = "CM Med trait vs micro.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
