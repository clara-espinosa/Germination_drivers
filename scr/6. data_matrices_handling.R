library(ggplot2);library(tidyverse);library(rstatix);library(lubridate);library (vegan)
# Data matrices handling
############################################# MEDITERRANEAN ###########################################################
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
unique(spatial_sp_Med$species)
setdiff(spatial_sp_Med$species, species$species)
# 2A- Cover/plot of species with germination traits ####
species %>%
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Med, by = "species")%>%
  filter (!species == "Solidago virgaurea")%>% # species with 0 germ
  filter (!species == "Teedalia conferta")%>% 
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover))%>%
  filter(rel_cover>79)->med_plot1
as.vector(med_plot1[,1])->med_plot65
# 15 plots with more than 80% coverage with sp traits (drivers) (considering raw cover, not to 100%)
# 65 plots with more than 80% coverage with sp traits (drivers) (considering relative cover, up to 100%)
# 64 plots with more than 80% coverage with sp traits (drivers + phenology)

# 2B- Number of species with germination traits per plot ####
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

# 3- temperature graphs of plots with 80% coverage with trait data NEED TO UPDATE ####
med_plot1%>%
  merge(read.csv("data/spatial-survey-temperatures-Med.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Med.csv"))-> spatial1 #%>%  #
#filter(site== "Rabinalto")-> spatial1

setdiff(med_plot1$plot, spatial1$plot)
unique(med_plot1$plot)
unique(spatial1$plot)

spatial1 %>% 
  group_by(Time = lubridate::floor_date(Time, "day"), site, plot) %>%
  summarise(Temperature = mean(Temperature)) -> spatial2

x11()
spatial1 %>%
  mutate(site = fct_relevel(site,"Rabinalto", "Cañada","Solana","Penouta")) %>%
  ggplot(aes(x=Time, y=Temperature, color = site)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = .05, alpha = 0.5) +
  geom_line(data = spatial2, linewidth = .75) +
  facet_wrap(~ plot,nrow = 10 ) + #   ,  labeller = plot
  labs(title = "Spatial survey (Jul 2021 - May 2022)", subtitle = "Plots with 80% cover with germination driver traits", 
       x= "Time (4-h recording inverval)", y= "Temperature (ºC)") +
  ggthemes::theme_tufte() +
  coord_cartesian(ylim = c(-5, 45)) +
  #scale_color_manual( values = c( "limegreen")) +
  scale_color_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4")) +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        legend.position = c(0.75,0.025), # 
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 12), #, face = "italic"
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 11, color = "black")) 

# DATA MATRICES ####
# header species data ####
read.csv("data/species.csv") %>%
  dplyr::select(species,code, family, community, habitat, germ_drivers)%>%
  filter(community == "Mediterranean")%>% # no germination
  filter(!species == "Teesdalia conferta") %>% # all germinated in cold stratification
  filter(!species=="Solidago virgaurea")-> spe_med
unique(spe_med$species)
# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
spe_med$sp <- make.cepnames(spe_med$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x plot matrix #####
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  filter(!species=="Solidago virgaurea")%>%
  filter(!species == "Teesdalia conferta") %>%
  merge(spe_med, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  spread(species, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  filter(!plot=="A00")%>% # plot with Microlog data NO iButton
  filter(!plot=="B00")%>%# plot with Microlog data NO iButton
  filter(!plot=="C00")%>%# plot with Microlog data NO iButton
  filter(!plot=="D00")%>%# plot with Microlog data NO iButton
  filter(!plot=="B13")%>% # iButtons could no be recovered
  filter(!plot=="B14")%>% # iButtons could no be recovered
  filter(!plot=="B15")%>% # iButtons could no be recovered
  filter(!plot=="C04")%>% # iButtons could no be recovered
  filter(!plot=="C07")%>% # iButtons could no be recovered
  filter(!plot=="D06")%>% # filter plots without vegetation
  filter(!plot=="D07")%>% # filter plots without vegetation
  merge(med_plot65, by= "plot")%>% # plots with more that 80% coverage with traits
  column_to_rownames(var="plot")-> plot_x_sp_M # replace Na with 0
sp_x_plot_M <- t(plot_x_sp_M)

# presence_absence sp x plot matrix (TO USE in CM)
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  filter(!species=="Solidago virgaurea")%>%
  filter(!species == "Teesdalia conferta") %>%
  merge(spe_med, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  filter(!plot=="A00")%>% # filter plots without iButtons data or without any vegetation
  filter(!plot=="B00")%>%
  filter(!plot=="C00")%>%
  filter(!plot=="D00")%>%
  filter(!plot=="B13")%>%
  filter(!plot=="B14")%>%
  filter(!plot=="B15")%>%
  filter(!plot=="C04")%>%
  filter(!plot=="C07")%>%
  merge(med_plot65, by= "plot")%>%
  column_to_rownames(var="plot")-> plot_x_sp_M_PA 

setdiff(med_plot65$plot, test1$plot)
# 2-  species x traits matrix ########
#### species x traits matrix with germination as log odd ratios 
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
  filter (!species == "Solidago virgaurea")%>%  #filter species with 0 germination traits
  filter(!species == "Teesdalia conferta") %>% # all germinated in pretreatment
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  mutate(treatment = as.factor(treatment))%>%
  merge(spe_med)-> germ_to_odds_M
unique(germ_to_odds_M$species)
str(germ_to_odds_M)

### Get GLM coefficients since the coefficients(estimate) are returned in log odds, 
#https://www.r-bloggers.com/2016/11/introducing-r-package-oddsratio/

glms <- function(x) {
  glm(cbind(finalgerm, viable - finalgerm) ~ treatment,  family = "binomial", data=x) -> m1
  broom::tidy(m1)
}

germ_to_odds_M %>%
  group_by(code, species) %>%
  do(glms(.)) %>%
  dplyr:: select( species, term, estimate) %>%
  rename(log_odds_ratio = estimate ) %>%
  mutate(term = as.factor(term))%>%
  mutate(term = fct_recode(term, odds_A_control = "(Intercept)", odds_B_dark = "treatmentB_alternate_dark",
                           odds_C_WP = "treatmentC_alternate_WP", odds_D_constant = "treatmentD_constant_light"))%>%
  dplyr::group_by(species)%>%
  tidyr::spread(term, log_odds_ratio)%>%
  merge(spe_med, by= c("species", "code"))-> germ_odds_M
unique(germ_odds_M$species)

# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  mutate(seed_mass= log (seed_mass))%>%
  right_join(spe_med, by= "species") -> seed_mass_M 
unique(seed_mass_M$species)
hist(seed_mass_M$seed_mass)

# plant height (all sp covered)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2))%>%
  mutate(plant_height=log(plant_height))%>%
  right_join(spe_med, by= "species")-> plant_height_M 
unique (plant_height_M$species)
hist(plant_height_M$plant_height)

#leaves traits (leaves traits measured x individual ) Jasione laevis missing
read.csv("data/leaves_traits.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species, N_leaves) %>%
  dplyr::summarise(leaf_area= round(mean(leaf_area),2),
                   LDMC = round(mean(LDMC),2),
                   SLA = round(mean(SLA),2))%>%
  group_by (species)%>%
  dplyr::summarise(leaf_area= leaf_area/N_leaves,
                   LDMC = LDMC/N_leaves,
                   SLA = SLA/N_leaves)%>%
  mutate(leaf_area= log (leaf_area),
          LDMC = LDMC,
          SLA = SLA)%>%
  right_join(spe_med, by= "species")-> leaves_traits_M 
unique (leaves_traits_M $species)
hist(leaves_traits_M $leaf_area)
hist(leaves_traits_M $LDMC)
hist(leaves_traits_M $SLA)

# join traits subset
germ_odds_M%>% # log odds ratios
  merge(seed_mass_M)%>% # log transformed
  merge(plant_height_M)%>% # log transformed
  merge(leaves_traits_M )%>% #only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height,leaf_area, LDMC, SLA,
                odds_B_dark, odds_C_WP, odds_D_constant)%>% # percentage data not necesary in theory
  column_to_rownames(var="species")->sp_x_trait_M

# 3-  plot x Environmental data ########
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
  dplyr::select(plot, elevation,bio1:GDD)%>% 
  filter(plot%in%med_plot65$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_env_M



################################# TEMPERATE ########################################################################
# 1-check species names to match####
# dataframe with species data from germination experiments
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers)-> species 

# dataframe with species data from spatial survay
read.csv("data/spatial-survey-species-Tem.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Tem
setdiff(spatial_sp_Tem$species, species$species)

# 2A- Cover/plot of species with germination traits ####
species %>%
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Tem, by = "species")%>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover))%>%
  filter(rel_cover>79)->tem_plot1
as.vector(tem_plot1[,1])->tem_plot48
# 1 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)
# 48 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)
# 40 plot with more than 80% coverage with sp traits (drivers + phenology) (considering raw cover, not to 100%)

# 2B- Number of species with germination traits per plot ####
species %>%
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Tem, by = "species")%>%
  #filter (!species == "Euphrasia salisburgensis")%>% # species with 0 germ across all treatments
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

# 3- temperature graphs of plots with 80% coverage with trait data NEED TO UPDATE ####
tem_plot1%>%
  merge(read.csv("data/spatial-survey-temperatures-Tem.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Tem.csv"))-> spatial1 

setdiff(tem_plot1$plot, spatial1$plot)
unique(tem_plot1$plot)
unique(spatial1$plot)

spatial1 %>% 
  group_by(Time = lubridate::floor_date(Time, "day"), site, plot) %>%
  summarise(Temperature = mean(Temperature)) -> spatial2

x11()
spatial1 %>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra")) %>%
  ggplot(aes(x=Time, y=Temperature, color = site)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = .05, alpha = 0.5) +
  geom_line(data = spatial2, linewidth = .75) +
  facet_wrap(~ plot,nrow = 10 ) + #   ,  labeller = plot
  labs(title = "Spatial survey (Oct 2018 - Aug 2019)", subtitle = "Plots with 80% cover with germination driver traits", 
       x= "Time (4-h recording inverval)", y= "Temperature (ºC)") +
  ggthemes::theme_tufte() +
  coord_cartesian(ylim = c(-5, 45)) +
  #scale_color_manual( values = c( "limegreen")) +
  scale_color_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4")) +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        legend.position = c(0.75,0.025), # 
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = 12), #, face = "italic"
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 11, color = "black")) 


# DATA MATRICES ####
# header species data  #######
read.csv("data/species.csv") %>%
  dplyr::select(species, code, family, community, habitat, germ_drivers)%>%
  filter(community == "Temperate")%>%
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Minuartia CF")-> spe_tem
# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
spe_tem$sp <- make.cepnames(spe_tem$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x plot matrix #########
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species))) %>%
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Minuartia CF")%>%
  merge(spe_tem, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  spread(species, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>%
  merge(tem_plot48, by= "plot")%>% # plots with more that 80% coverage with traits
  column_to_rownames(var="plot")-> plot_x_sp_T # replace Na with 0
sp_x_plot_T <- t(plot_x_sp_T)
setdiff(spe_tem$species, test1$species)

# presence_absence sp x plot matrix (to use in CM)
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Minuartia CF")%>%
  merge(spe_tem, by="species")%>%  # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% 
  merge(tem_plot48, by= "plot")%>%
  column_to_rownames(var="plot")-> plot_x_sp_T_PA 
# 2-  species x traits matrix ###########
#### species x traits matrix with germination  as odd ratios 
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
  filter (!species == "Euphrasia salisburgensis")%>%  #filter species with 0 germination traits?
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  merge(spe_tem)-> germ_to_odds_T

### Get GLM coefficients since the coefficients(estimate) are returned in log odds,
#https://www.r-bloggers.com/2016/11/introducing-r-package-oddsratio/
glms <- function(x) {
  glm(cbind(finalgerm, viable - finalgerm) ~ treatment,  family = "binomial", data= x) -> m1
  broom::tidy(m1)
}

germ_to_odds_T %>%
  group_by(code, species) %>%
  do(glms(.)) %>%
  dplyr:: select(code, species, term, estimate) %>%
  rename(log_odds_ratio = estimate ) %>%
  mutate(term = as.factor(term))%>%
  mutate(term = fct_recode(term, odds_A_control = "(Intercept)", odds_B_dark = "treatmentB_alternate_dark",
                           odds_C_WP = "treatmentC_alternate_WP", odds_D_constant = "treatmentD_constant_light"))%>%
  dplyr::group_by(species)%>%
  tidyr::spread(term, log_odds_ratio)%>%
  merge(spe_tem, by= c("species", "code"))-> germ_odds_T

# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  mutate(seed_mass=log(seed_mass))%>%
  merge(spe_tem, by= "species") -> seed_mass_T 
hist(seed_mass_T $seed_mass) 

# plant and floral height (only Gentianella campestris missing but is excluded due to no germination across Treat)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2))%>%
  mutate(plant_height=log(plant_height))%>%
  right_join(spe_tem, by= "species")-> plant_height_T

hist(plant_height_T$plant_height) 
 
#leaves traits (leaves traits measured x individual ) all species covered
read.csv("data/leaves_traits.csv", sep= ",") %>%
  filter(community == "Temperate")%>%
  group_by(species, N_leaves) %>%
  dplyr::summarise(leaf_area= round(mean(leaf_area),2),
                   LDMC = round(mean(LDMC),2),
                   SLA = round(mean(SLA),2))%>%
  group_by(species)%>%
  dplyr::summarise(leaf_area= leaf_area/N_leaves,
                   LDMC = LDMC/N_leaves,
                   SLA = SLA/N_leaves)%>%
  mutate(leaf_area= log (leaf_area),
         LDMC = LDMC,
         SLA = SLA)%>%
  right_join(spe_tem, by= "species")-> leaves_traits_T 
unique (leaves_traits_T $species)
hist(leaves_traits_T$leaf_area)
hist(leaves_traits_T$LDMC)
hist(leaves_traits_T$SLA)


# join traits subset
germ_odds_T%>% # log odd ratios
  merge(seed_mass_T)%>% #log seed mass
  merge(plant_height_T)%>% # log plant height
  merge(leaves_traits_T)%>% # only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height, 
                leaf_area, LDMC, SLA, 
                odds_B_dark, odds_C_WP, odds_D_constant)%>%
  column_to_rownames(var="species")->sp_x_trait_T

# 3-  plot x Environmental data ##########
read.csv("data/spatial-survey-temperatures-Tem.csv") %>%
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
  merge(read.csv("data/spatial-survey-header-Tem.csv")) %>% 
  dplyr::select(plot, elevation,bio1:GDD)%>%
  #dplyr::select(plot, elevation, FDD, GDD)%>%
  filter(plot%in%tem_plot48$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_env_T
