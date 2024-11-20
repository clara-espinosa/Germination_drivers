library(ggplot2);library(tidyverse);library(rstatix);library(lubridate);library (vegan)
# Data matrices handling
############################################# MEDITERRANEAN ###########################################################
# 1-check species names to match####
# dataframe with species data from germination experiments
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers) %>%
  filter (community == "Mediterranean")%>%
  filter(!species=="Solidago virgaurea")%>% # 0 germ across all treatments
  filter(!species=="Avenella flexuosa")%>% # 0 germ across all treatments
  filter(!species=="Cerastium ramosissimum")%>% # 0 germ across all treatments
  filter(!species=="Phalacrocarpum oppositifolium")%>% # 0 germ across all treatments
  filter(!species == "Teesdalia conferta") %>% # all germinated during cold stratification
  filter(!species == "Jasione laevis")%>% # without leaves traits
  as.data.frame()-> sp_med 
unique(sp_med$species) # 18 species with full traits
# dataframe with species data from spatial survay
read.csv("data/spatial-survey-species-Med.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Med
unique(spatial_sp_Med$species)
setdiff(species$species,spatial_sp_Med$species) 
# 2- Cover/plot of species with germination traits ####
sp_med %>%
  merge(spatial_sp_Med, by = "species")%>%
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover),
            N_sp = length(unique(species)))%>%
  filter(rel_cover>79)%>%
  filter(N_sp>2)%>%
  dplyr:: select(plot, N_sp)->med_plot

# 53 plots with more than 80% coverage with sp traits (drivers) (considering relative cover, up to 100%)

# 3- temperature graphs of plots with 80% coverage with trait data ####
med_plot%>%# 53 plots with >80% coverage with trait and environmental data
  merge(read.csv("data/spatial-survey-temperatures-Med.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Med.csv"))%>%
  filter(!plot=="A00")%>% # plot with Microlog data NO iButton
  filter(!plot=="B00")%>%# plot with Microlog data NO iButton
  filter(!plot=="C00")%>%# plot with Microlog data NO iButton
  filter(!plot=="D00")%>%# plot with Microlog data NO iButton
  filter(!plot=="B13")%>% # iButtons could no be recovered
  filter(!plot=="B14")%>% # iButtons could no be recovered
  filter(!plot=="B15")%>% # iButtons could no be recovered
  filter(!plot=="C07")-> spatial_env_med 

setdiff(med_plot$plot, spatial_env_med$plot)
setdiff(spatial_env_med$plot,med_plot$plot)
unique(med_plot$plot)
unique(spatial_env_med$plot)

spatial_env_med %>% 
  group_by(Time = lubridate::floor_date(Time, "day"), site, plot) %>%
  summarise(Temperature = mean(Temperature)) -> spatial_mean_med

x11()
spatial_env_med %>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  ggplot(aes(x=Time, y=Temperature, color = site)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = .05, alpha = 0.5) +
  geom_line(data = spatial_mean_med, linewidth = .75) +
  facet_wrap(~ plot,nrow = 10 ) + #   ,  labeller = plot
  labs(title = "Spatial survey (July 2021 - May 2022)", subtitle = "Plots with 80% relative cover with traits", 
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
        axis.text.y = element_text(size = 11, color = "black"))->tem_plots_Med;tem_plots_Med

ggsave(filename = "Temperatures plots with traits Med.png", plot =tem_plots_Med , path = "results/Supplementary/", 
       device = "png", dpi = 600)

# DATA MATRICES ####
# header species data ####
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers) %>%
  filter (community == "Mediterranean")%>%
  filter(!species=="Solidago virgaurea")%>% # 0 germ across all treatments
  filter(!species == "Teesdalia conferta") %>% # all germinated during cold stratification
  filter(!species == "Jasione laevis")%>% # without leaves traits
  as.data.frame()-> sp_med 
unique(sp_med$species) # 21 species with full traits
# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
sp_med$sp <- make.cepnames(sp_med$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x plot matrix #####
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  merge(sp_med, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  spread(species, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>% #-> med_plot
  filter(plot%in%spatial_env_med$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_sp_M # replace Na with 0
dim(plot_x_sp_M)
sp_x_plot_M <- t(plot_x_sp_M)
dim(plot_x_sp_M)
unique(sp_med$species)

# presence_absence sp x plot matrix (to use in CM) with restrictive dataset
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  merge(sp_med, by="species")%>%  # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% 
  filter(plot%in%spatial_env_med$plot)%>%
  column_to_rownames(var="plot")-> plot_x_sp_M_PA

sp_x_plot_M <- t(plot_x_sp_M)
dim(plot_x_sp_M)
dim(plot_x_sp_M_PA)

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
  filter(!species=="Avenella flexuosa")%>% # 0 germ across all treatments
  filter(!species=="Cerastium ramosissimum")%>% # 0 germ across all treatments
  filter(!species=="Phalacrocarpum oppositifolium")%>% # 0 germ across all treatments
  filter(!species == "Teesdalia conferta") %>% # all germinated in pretreatment
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  mutate(treatment = as.factor(treatment))%>%
  merge(sp_med)-> germ_to_odds_M
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
  merge(sp_med, by= c("species", "code"))-> germ_odds_M
unique(germ_odds_M$species)

# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  mutate(seed_mass= log (seed_mass))%>% # almost normal ??
  right_join(sp_med, by= "species") -> seed_mass_M 
unique(seed_mass_M$species)
hist(seed_mass_M$seed_mass)

# plant height (all sp covered)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Mediterranean")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2))%>%
  mutate(plant_height=log(plant_height))%>% # looks almost normal
  right_join(sp_med, by= "species")-> plant_height_M 
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
  right_join(sp_med, by= "species")-> leaves_traits_M 
unique (leaves_traits_M $species)
hist(leaves_traits_M$leaf_area) # all look normal
hist(leaves_traits_M$LDMC)
hist(leaves_traits_M$SLA)# 

# join traits subset
germ_odds_M%>% # log odds ratios
  merge(seed_mass_M)%>% # log transformed
  merge(plant_height_M)%>% # log transformed
  merge(leaves_traits_M )%>% #only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height,leaf_area, LDMC, SLA,
                odds_B_dark, odds_C_WP, odds_D_constant)%>% # percentage data not necesary in theory
  column_to_rownames(var="species")->sp_x_trait_M
dim(sp_x_trait_M)

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
  filter(plot%in%spatial_env_med$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_env_M
dim(plot_x_env_M)
plot_x_env_M%>%
  rownames_to_column(var = "plot")->plot_x_env_M2 
plot_x_env_M2 %>%
  get_summary_stats()

################################# TEMPERATE ########################################################################
# 1-check species names to match####
# dataframe with species data from germination experiments
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers)%>%
  filter (community == "Temperate") %>%
  filter (!species == "Euphrasia salisburgensis")%>%  #0 germ across all treatments
  filter (!species == "Gentiana verna")%>% #0 germ across all treatments
  filter (!species == "Gentianella campestris")%>% #0 germ across all treatments
  filter (!species == "Kobresia myosuroides")%>% #0 germ across all treatments
  filter (!species == "Salix breviserrata")%>% #0 germ across all treatments
  filter (!species == "Sedum album")%>% #0 germ across all treatments
  filter (!species == "Sedum atratum")%>% #0 germ across all treatments
  filter (! species == "Veronica nummularia") %>% #0 germ across all treatments
  filter (!species == "Minuartia CF")%>% #missing traits 
  as.data.frame()-> sp_tem 
unique(sp_tem$species) # 25 species with full traits
# dataframe with species data from spatial survay
read.csv("data/spatial-survey-species-Tem.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Tem
setdiff(spatial_sp_Tem$species, species$species)

# 2- Cover/plot of species with germination traits ####
sp_tem %>%
  merge(spatial_sp_Tem, by = "species")%>%
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover),
            N_sp = length(unique(species)))%>%
  filter(rel_cover>79)%>%
  filter(N_sp>2)%>%
  dplyr:: select(plot, N_sp)->tem_plot

# 47 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)

# 3- temperature graphs of plots with 80% coverage with trait data NEED TO UPDATE ####
tem_plot%>%
  merge(read.csv("data/spatial-survey-temperatures-Tem.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"))-> spatial_env_tem 

setdiff(tem_plot$plot, spatial_env_tem$plot)
unique(tem_plot$plot)
unique(spatial_env_tem$plot)

spatial_env_tem %>% 
  group_by(Time = lubridate::floor_date(Time, "day"), site, plot) %>%
  summarise(Temperature = mean(Temperature)) -> spatial_mean_tem

x11()
spatial_env_tem  %>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra")) %>%
  ggplot(aes(x=Time, y=Temperature, color = site)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = .05, alpha = 0.5) +
  geom_line(data = spatial_mean_tem, linewidth = .75) +
  facet_wrap(~ plot,nrow = 10 ) + #   ,  labeller = plot
  labs(title = "Spatial survey (Oct 2018 - Aug 2019)", subtitle = "Plots with 80% relative cover with traits", 
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
        axis.text.y = element_text(size = 11, color = "black")) ->tem_plots_Tem;tem_plots_Tem

ggsave(filename = "Temperatures plots with traits Tem.png", plot =tem_plots_Tem , path = "results/Supplementary/", 
       device = "png", dpi = 600)

# DATA MATRICES ####
# header species data  #######
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers)%>%
  filter (community == "Temperate") %>%
  filter (!species == "Euphrasia salisburgensis")%>%  #0 germ across all treatments
  filter (!species == "Gentiana verna")%>% #0 germ across all treatments
  filter (!species == "Gentianella campestris")%>% #0 germ across all treatments
  filter (!species == "Kobresia myosuroides")%>% #0 germ across all treatments
  filter (!species == "Salix breviserrata")%>% #0 germ across all treatments
  filter (!species == "Sedum album")%>% #0 germ across all treatments
  filter (!species == "Sedum atratum")%>% #0 germ across all treatments
  filter (!species == "Veronica nummularia")%>% #0 germ across all treatments
  filter (!species == "Minuartia CF")%>% #missing traits 
  as.data.frame()-> sp_tem 
unique(sp_tem$species) # 26 species with full traits
# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
sp_tem$sp <- make.cepnames(sp_tem$species, seconditem=FALSE) #works if Sp are in columns

# 1-  species x plot matrix #########
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species))) %>%
  merge(sp_tem, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  spread(species, rel_cover)%>% # species in columns
  replace(is.na(.), 0)%>%
  filter(plot%in%spatial_env_tem$plot)%>% # plots with more that 80% coverage with traits
  column_to_rownames(var="plot")-> plot_x_sp_T # replace Na with 0
sp_x_plot_T <- t(plot_x_sp_T)

# presence_absence sp x plot matrix (to use in CM) with restrictive dataset
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  merge(sp_tem, by="species")%>%  # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% 
  filter(plot%in%spatial_env_tem$plot)%>%
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
  filter (!species == "Veronica nummularia")%>%
  filter (!species == "Minuartia CF")%>%
  spread(treatment, finalgerm)%>%
  mutate (B_alternate_dark = B_alternate_dark - E_cold_stratification)  %>% # substract cold stratification
  mutate (B_alternate_dark = ifelse(B_alternate_dark < 0, 0, B_alternate_dark)) %>% # convert negative values to 0 germination
  gather ("treatment", "finalgerm", A_alternate_light:E_cold_stratification) %>%
  filter(!treatment=="E_cold_stratification")%>%
  merge(viables_sp, by= c("code", "species", "treatment"))%>%
  merge(sp_tem)-> germ_to_odds_T

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
  merge(sp_tem, by= c("species", "code"))-> germ_odds_T

# seed mass trait (all sp covered)
read.csv("data/seed_mass.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(seed_mass= round(mean(weight_50_mg),2))%>%
  mutate(seed_mass=log(seed_mass))%>%
  merge(sp_tem, by= "species") -> seed_mass_T 
hist(seed_mass_T $seed_mass) # looks normal

# plant and floral height (only Gentianella campestris missing but is excluded due to no germination across Treat)
read.csv("data/plant_height.csv", sep = ",") %>%
  filter(community == "Temperate")%>%
  group_by(species) %>%
  dplyr::summarize(plant_height= round(mean(plant_height),2))%>%
  mutate(plant_height=log(plant_height))%>%
  right_join(sp_tem, by= "species")-> plant_height_T

hist(plant_height_T$plant_height) # looks normal distributed
 
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
  right_join(sp_tem, by= "species")-> leaves_traits_T 
unique (leaves_traits_T $species)
hist(leaves_traits_T$leaf_area) # look almost normal
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
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep= ";")) %>% 
  dplyr::select(plot, elevation,bio1:GDD)%>%
  filter(plot%in%spatial_env_tem$plot)%>%
  column_to_rownames(var="plot")-> plot_x_env_T

plot_x_env_T%>%
  get_summary_stats()

plot_x_env_T%>%
  rownames_to_column(var = "plot")->plot_x_env_T2 
plot_x_env_T2 %>%
  get_summary_stats()

