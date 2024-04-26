library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych)
Sys.setlocale("LC_TIME", "English")

# preliminary steps #
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
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Med, by = "species")%>%
  filter (!species == "Solidago virgaurea")%>% # species with 0 germ
  group_by(plot) %>%
  #summarise(cover = sum(cover))%>%
  #filter(cover>79)%>%
  summarise(rel_cover = sum(rel_cover))%>%
  filter(rel_cover>79)->med_plot1
# 15 plots with more than 80% coverage with sp traits (drivers) (considering raw cover, not to 100%)
# 65 plots with more than 80% coverage with sp traits (drivers) (considering relative cover, up to 100%)
# 64 plots with more than 80% coverage with sp traits (drivers + phenology)

#TEMPERATE
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
# 1 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)
# 48 plot with more than 80% coverage with sp traits (drivers)(considering raw cover, not to 100%)
# 40 plot with more than 80% coverage with sp traits (drivers + phenology) (considering raw cover, not to 100%)

# 2B- Number of species with germination traits per plot ####
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
# 3- temperature graphs of plots with 80% coverage with trait data ####
# MEDITERRANEAN
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

# Temperate
tem_plot1%>%
  merge(read.csv("data/spatial-survey-temperatures-Tem.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Tem.csv"))-> spatial1 #%>%  #
#filter(site== "Rabinalto")-> spatial1

setdiff(tem_plot1$plot, spatial1$plot)
unique(tem_plot1$plot)
unique(spatial1$plot)

spatial1 %>% 
  group_by(Time = lubridate::floor_date(Time, "day"), site, plot) %>%
  summarise(Temperature = mean(Temperature)) -> spatial2

unique(spatial2$site)
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


# CWM test  Following chapter 5 R TRait-Based ecology by Lars####
# 5.1 Data matrices ####
# base matrices from species level analysis but with reformatting
# in sp x plot (species in rows and plots in columns)
# in all files names of species and plots are not part of the matrix, used as row and columns names
# mediterranean #
sp_x_com_M <- t(sp_x_com_M) # abundance already as relative cover
sp_x_trait_M
com_x_env_M

# temperate #
sp_x_com_T <- t(sp_x_com_T)# abundance already as relative cover
sp_x_trait_T
com_x_env_T

## 5.2 Check normality of the quantitative traits 
### all traits are right skewed, ty a log transformation look better but still not normal
sp_x_trait_M%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

sp_x_trait_T%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

## 5.3 Calculation of CWM ####
# Mediterranean
CWM_M<-functcomp(sp_x_trait_M, t(sp_x_com_M), CWM.type = "all") # CWM.type to consider binary or categorical traits
# Let’s see if CWM changes along the climatic gradient
# explore visually a bit the results before running any analysis,
com_x_env_M%>%
  rownames_to_column(var = "plot")->com_x_env_M2 

CWM_M%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  merge(com_x_env_M2)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(trait = factor(trait))%>%
  #mutate(trait = recode (trait, "Immediate" = "Fresh", "After_ripening" = "After ripened"))%>%
  mutate(trait = fct_relevel(trait, "A_alternate_light","B_alternate_dark","C_alternate_WP","D_constant_light", "E_cold_stratification",
                             "odds_A_control", "odds_B_dark","odds_C_WP","odds_D_constant","odds_E_cold", 
                             "plant_height", "floral_height", "seed_mass", "seed_production"))%>%
  ggplot(aes(y=value, x= elevation, fill = site), color= "black")+
  geom_point(shape = 21, size =3)+
  geom_smooth (aes(y=value, x= elevation),method = "lm", se=FALSE, color = "black", inherit.aes=F)+
  facet_wrap(~trait, nrow = 3, scales = "free")+
  labs(title="Community Weighted Means in Mediterranean system")+
  theme_classic(base_size = 12)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")
# GDD, FDD, Snow, Elevation
# Temperate
CWM_T<-functcomp(sp_x_trait_T, t(sp_x_com_T), CWM.type = "all") # CWM.type to consider binary or categorical traits
# Let’s see if CWM changes along the climatic gradient
# explore visually a bit the results before running any analysis,
com_x_env_T%>%
  rownames_to_column(var = "plot")->com_x_env_T2 

CWM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  merge(com_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(trait = factor(trait))%>%
  #mutate(trait = recode (trait, "Immediate" = "Fresh", "After_ripening" = "After ripened"))%>%
  mutate(trait = fct_relevel(trait, "A_alternate_light","B_alternate_dark","C_alternate_WP","D_constant_light", "E_cold_stratification",
                             "odds_A_control", "odds_B_dark","odds_C_WP","odds_D_constant","odds_E_cold", 
                             "plant_height", "floral_height", "seed_mass", "seed_production"))%>%
  ggplot(aes(y=value, x= Snw, fill = site), color= "black")+
  geom_point(shape = 21, size =3)+
  geom_smooth (aes(y=value, x= Snw),method = "lm", se=FALSE, color = "black", inherit.aes=F)+
  facet_wrap(~trait, nrow = 3, scales = "free")+
  labs(title="Community Weighted Means in Temperate system")+
  theme_classic(base_size = 12)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")
# GDD, FDD, Snow, Elevation
# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(CWM_T$E_cold_stratification~elevation+ FDD + GDD + Snw, data=com_x_env_T))
# see summary in results
# we can also use the multivariate analyses for the representation, for example RDA

 rda_med_all<-rda(CWM_M~elevation , data= com_x_env_M)
 plot(rda_med_all, type = "n", scaling = "sites")
 text(rda_med_all, dis = "cn", scaling = "sites")
 text(rda_med_all, dis = "sp", scaling = "sites", col = "red")
 