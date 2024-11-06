library(ggplot2);library(tidyverse);library(rstatix)
library(lubridate);library (vegan);library(dplyr);library(vegan)
library(lubridate);library(FD); library(psych);library(picante)
library (gawdis)

####### MEDITERRANEAN ######
# Number of species with germination traits per plot (extended data set) ####
sp_med %>%
  merge(spatial_sp_Med, by = "species")%>%
  select(species, plot, N_sp)%>%
  group_by(plot) %>%
  summarise(N_sp = first(N_sp), 
            sp_traits = length(unique(species)), 
            sp_covered = sp_traits/N_sp)%>%
  filter(sp_covered>0.49) %>%
  filter(sp_traits>1)->med_plot_extended
#arrange(sp_covered)%>%
print(n=84)
# in 44 plots more than 75% sp covered with drivers traits (focused on presence)
# in 70 plots more than 49% sp covered with drivers traits and more than 2 sp per plot
# Temperature graphs of plot with more than 50% species with traits and at least 2 species ####
med_plot_extended%>%# 70 plots 
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
  filter(!plot=="C07")-> spatial_env_med_extended

setdiff(med_plot_extended$plot, spatial_env_med_extended$plot)
setdiff(spatial_env_med_extended$plot,med_plot_extended$plot)
unique(med_plot_extended$plot)
unique(spatial_env_med_extended$plot)

spatial_env_med_extended %>% 
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
# Presence_absence sp x plot matrix (TO USE in CM) with extended data set ####
read.csv("data/spatial-survey-species-Med.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  merge(sp_med, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% 
  filter(plot%in%spatial_env_med_extended$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_sp_M_PA_ext 

# Plot x Environmental data extended dataset for CM ####
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
  merge(disturbance_M, by="plot")%>%
  filter(plot%in%spatial_env_med_extended$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_env_M_ext

plot_x_env_M_ext%>%
  rownames_to_column(var = "plot")->plot_x_env_M2_ext 
plot_x_env_M2_ext %>%
  get_summary_stats()
# 5.1 Data matrices ####
# extended dataset for CM
sp_x_trait_M
plot_x_env_M_ext
plot_x_sp_M_PA_ext
# 5.2 Calculation of CM extended dataset ####
sp_x_trait_M<- as.matrix(sp_x_trait_M)
plot_x_sp_M_PA_ext<- as.matrix(plot_x_sp_M_PA_ext)

CM_M_ext<-functcomp(sp_x_trait_M, plot_x_sp_M_PA_ext, CWM.type = "all") # CWM.type to consider binary or categorical traits

# first let's check CM normality (all look normal distribution)
x11()
CM_M%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# all traits look quite normally distributed!!
plot_x_env_M%>%
  rownames_to_column(var = "plot")->plot_x_env_M2 

CM_M_ext%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_M2_ext)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CM_microclima_M_ext

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model 
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cm <- function(x) {
  glm(CM ~ scale(elevation) + scale(FDD) + scale(GDD) + scale(Snw), data = x) -> m1
  broom::tidy(m1)
}

CM_M_ext%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_M2_ext, by= "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cm(.)) %>%
  write.csv("results/CM ext vs micro MED.csv")

# double facetting with significances
ann.sig.CM.M <- data.frame (read.csv("results/sig_CM_M.csv", sep = ";"))
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                 "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")
str(CM_microclima_M_ext)
x11()
CM_microclima_M_ext%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(factor(trait, levels = c("odds_B_dark","odds_C_WP","odds_D_constant",
                                      "seed_mass","plant_height", 
                                      "leaf_area", "LDMC", "SLA"))~micro_variable, scales = "free", 
             labeller = as_labeller(trait_names))+
  geom_text(data=ann.sig.CM.M, label =ann.sig.CM.M$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Community Means in Mediterranean system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> plot_CM_M;plot_CM_M
ggsave(plot_CM_M, file = "CM Med trait vs micro.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
#width = 320, height = 260, units = "mm", 

# see summary in results
# we can also use the multivariate analyses for the representation, for example RDA

rda_med_all<-rda(CM_M_ext~elevation+ FDD + GDD + Snw , data= plot_x_env_M_ext)
plot(rda_med_all, type = "n", scaling = "sites")
text(rda_med_all, dis = "cn", scaling = "sites")
text(rda_med_all, dis = "sp", scaling = "sites", col = "red")

rda_med_0<-rda(CM_M_ext~1, data= plot_x_env_M_ext)
rda_med_all<-rda(CM_M_ext~elevation+ FDD + GDD + Snw , data= plot_x_env_M_ext)
ordistep(rda_med_0, scope=formula(rda_med_all), direction = "forward")

RsquareAdj (rda_med_all)$adj.r.squared # 0.06
# final model only with elevation


# 5.3 GLMs CM x micro ####
read.csv("data/spatial-survey-header-Med.csv")%>%
  filter(plot%in%spatial_env_med_extended$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_M_ext # matrix with coordenates of Mediterranean plots

glms.auto.M.ext <- function(x) {
  glm1 <- glm(CM_ext ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_M_ext , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(CM_ext ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}

CM_M_ext%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, CM_ext, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_M2_ext) %>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, CM_ext,
                Snw, FDD, GDD, elevation)%>%
  group_by(trait)%>%
  do(glms.auto.M.ext(.))%>%
  write.csv("results/GLMs auto Med ext.csv")

####### TEMPERATE ######
# Number of species with germination traits per plot (extended data set) ####
sp_tem %>%
  merge(spatial_sp_Tem, by = "species")%>%
  select(species, plot, N_sp)%>%
  group_by(plot) %>%
  summarise(N_sp = first(N_sp), 
            sp_traits = length(unique(species)), 
            sp_covered = sp_traits/N_sp)%>%
  filter(sp_covered>0.49) %>%
  filter(sp_traits>1)-> tem_plot_extended
arrange(sp_covered)%>%
  print(n=84)
# only 5 plots more than 75% sp covered with drivers traits (focused on presence)
# 68 plots more than 49% of species covered with drivers traits
# Temperature graphs of plot with more than 50% species with traits and at least 2 species ####
tem_plot_extended%>%
  merge(read.csv("data/spatial-survey-temperatures-Tem.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Tem.csv"))-> spatial_env_tem_extended 

setdiff(tem_plot_extended$plot, spatial_env_tem_extended$plot)
unique(tem_plot_extended$plot)
unique(spatial_env_tem_extended$plot)

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

# Presence_absence sp x plot matrix (TO USE in CM) with extended data set ####
read.csv("data/spatial-survey-species-Tem.csv") %>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = round((cover/total_cover)*100,2), # relative cover according to max vegetation cover
         N_sp = length(unique(species)))%>%
  merge(sp_tem, by="species")%>% # get species names shortened
  dplyr::select(plot, species, rel_cover)%>%
  group_by(species, plot)%>%
  mutate(presence= ifelse(rel_cover>0,1,0))%>%
  dplyr::select(plot, species,presence)%>%
  spread(species, presence)%>% # species in columns
  replace(is.na(.), 0)%>% 
  filter(plot%in%spatial_env_tem_extended$plot)%>% 
  column_to_rownames(var="plot")-> plot_x_sp_T_PA_ext 
# Plot x Environmental dataextended dataset for CM ####
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
  merge(disturbance_T, by= "plot")%>%
  filter(plot%in%spatial_env_tem_extended$plot)%>%
  column_to_rownames(var="plot")-> plot_x_env_T_ext

plot_x_env_T_ext%>%
  get_summary_stats()

plot_x_env_T_ext%>%
  rownames_to_column(var = "plot")->plot_x_env_T2_ext
plot_x_env_T2_ext %>%
  get_summary_stats()

# 5.1 Data matrices ####
# extended dataset for CM
sp_x_trait_T
plot_x_env_T_ext
plot_x_sp_T_PA_ext

# 5.2 Calculation of CM  extended dataset ####
sp_x_trait_T<- as.matrix(sp_x_trait_T)
plot_x_sp_T_PA_ext<- as.matrix(plot_x_sp_T_PA_ext)

CM_T_ext<-functcomp(sp_x_trait_T, plot_x_sp_T_PA_ext, CWM.type = "all") # CWM.type to consider binary or categorical traits

# first let's check CWM normality (all look normal distribution)
CM_T_ext%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# all look pretty normal distributed

plot_x_env_T_ext%>%
  rownames_to_column(var = "plot")->plot_x_env_T2_ext

CM_T_ext%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_T2_ext)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep= ";"), by = c("plot", "elevation"))%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra")) %>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait,"odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass", "plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CM_microclima_T_ext

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cm <- function(x) {
  glm(CM ~ scale(elevation) + scale(FDD) + scale(GDD)+ scale(Snw), data = x) -> m1
  broom::tidy(m1)
}

CM_T_ext%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_T2_ext, by= "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cm(.)) %>%
  write.csv("results/CM ext vs micro TEM.csv")
# double facetting with significances
ann.sig.CM.T <- data.frame (read.csv("results/sig_CM_T.csv", sep = ";"))
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                 "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")

CM_microclima_T_ext%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(factor(trait, levels = c("odds_B_dark","odds_C_WP","odds_D_constant",
                                      "seed_mass","plant_height", 
                                      "leaf_area", "LDMC", "SLA"))~micro_variable, scales = "free", 
             labeller = as_labeller(trait_names))+
  geom_text(data=ann.sig.CM.T, label =ann.sig.CM.T$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Community Means in Temperate system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> plot_CM_T;plot_CM_T
ggsave(plot_CM_T, file = "CM Tem trait vs micro.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
#width = 320, height = 260, units = "mm", 


# we can also use the multivariate analyses for the representation, for example RDA

rda_tem_all<-rda(CM_T_ext~elevation+ FDD + GDD + Snw , data= plot_x_env_T_ext)
plot(rda_tem_all, type = "n", scaling = "sites")
text(rda_tem_all, dis = "cn", scaling = "sites")
text(rda_tem_all, dis = "sp", scaling = "sites", col = "red")

rda_tem_0<-rda(CM_T_ext~1, data= plot_x_env_T_ext)
rda_tem_all<-rda(CM_T_ext~elevation+ FDD + GDD + Snw , data= plot_x_env_T_ext)
ordistep(rda_tem_0, scope=formula(rda_tem_all), direction = "forward")

RsquareAdj (rda_tem_all)$adj.r.squared # 0.26
# final model with GDD and elevation significant and thus included!!
# 5.3 GLMs CM x micro ####
read.csv("data/spatial-survey-header-Tem.csv", sep= ";")%>%
  filter(plot%in%spatial_env_tem_extended$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_T_ext# matrix with coordenates of Temperate plots

glms.auto.T.ext <- function(x) {
  glm1 <- glm(CM_ext ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_T_ext , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(CM_ext ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}

CM_T_ext%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, CM_ext, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_T2_ext) %>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, CM_ext,
                Snw, FDD, GDD, elevation)%>%
  group_by(trait)%>%
  do(glms.auto.T.ext(.))%>%
  write.csv("results/GLMs auto Tem ext.csv")
