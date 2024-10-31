library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis); library(glmm); library(spdep) # para crear la autocovariate
library(glmmTMB); library (DHARMa) # paquetes para los GLMs

 ####################### MEDITERRANEAN ##########################################
# CM restricted dataset (from script 6)
read.csv("data/spatial-survey-header-Med.csv")%>%
  filter(plot%in%spatial_env_med$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_M # matrix with coordenates of Mediterranean plots

CM_M%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_M2, by= "plot")%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude, 
                elevation, Snw, FDD, GDD, 
                seed_mass:odds_D_constant)-> CM_M_glm

# process following Borja's script
# 1. model including latitude and longitude and the rest of env variables all scaled!!
glm1 <- glm(plant_height ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw), data = CM_M_glm)
# test with glmm (not same results)
# 2. Create autocovariate from "spdep" package from model results (step 4)
autocov <- autocov_dist(glm1$residuals, longlat_M , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE) #
# 3. add autocovariate to dataframe
CM_M_glm$auto <- autocov # add autocovariate variable to the dataset
# 3b. remove extreme values (only if error appears)
summary(CM_M_glm$auto) # para quitar valores extremos (daba error)
CM_M_glm <- subset(CM_M_glm, auto < Inf) 
CM_M_glm <- subset(CM_M_glm, auto > -500)

# 4. FINAL GLM- new model including the variable "auto" for accounting spatial autocorrelation
glm1.auto<- glm(plant_height ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = CM_M_glm)
summary(glm1.auto)
plot(glm1.auto)


### attemps to speed up the process #####
# try to combine data sets an run a loop group datasets with same plot subsets and coordinates
# meaning CM and CWM restricted dataset (separate per each community)
# CM_ext extended dataset apart
# FD apart due to add disturbance as explanatory variable

###### CM and CWM MEDITERRANEAN ####
read.csv("data/spatial-survey-header-Med.csv")%>%
  filter(plot%in%spatial_env_med$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_M # matrix with coordenates of Mediterranean plots

CM_M%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  merge(CWM_M%>%
    rownames_to_column(var = "plot")%>%
    gather(trait, CWM, seed_mass:odds_D_constant))%>%
  full_join(plot_x_env_M2) %>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, CM, CWM, 
                Snw, FDD, GDD, elevation)%>%
  gather(community_metric, value, CM:CWM)%>%
  group_by(trait, community_metric)%>%
  do(glms.auto.M(.))%>%
  write.csv("results/GLMs auto Med.csv")
  

glms.auto.M <- function(x) {
  glm1 <- glm(value ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_M , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(value ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}

###### CM_ext MEDITERRANEAN ####
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

## FD in Mediterranan ####
read.csv("data/spatial-survey-header-Med.csv")%>%
  filter(plot%in%spatial_env_med$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_M

glms.FD.M <- function(x) {
  glm1 <- glm(MPD ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw) + scale(disturbance), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_M , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(MPD ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(disturbance) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}
MPD.M(dark.dist.M, plot_x_sp_M)%>% # check function in script 7. 
  rbind(MPD.M(WP.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(Tconstant.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(germtraits.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(seedmass.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(plantheight.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(leafarea.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(LDMC.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(SLA.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(Wplanttraits.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(alltraits.dist.M, plot_x_sp_M))%>%
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.M", "" , trait))%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, data_type, MPD, 
                elevation,FDD, GDD, Snw,  disturbance)%>%
  group_by(trait, data_type)%>%
  do(glms.FD.M (.))%>%
  write.csv("results/GLMs auto FD Med.csv")
  
###### CM and CWM Temperate ####
read.csv("data/spatial-survey-header-Tem.csv", sep= ";")%>%
  filter(plot%in%spatial_env_tem$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_T# matrix with coordenates of Temperate plots

CM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  merge(CWM_T%>%
          rownames_to_column(var = "plot")%>%
          gather(trait, CWM, seed_mass:odds_D_constant))%>%
  full_join(plot_x_env_T2) %>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, CM, CWM, 
                Snw, FDD, GDD, elevation)%>%
  gather(community_metric, value, CM:CWM)%>%
  group_by(trait, community_metric)%>%
  do(glms.auto.T(.))%>%
  write.csv("results/GLMs auto Tem.csv")


glms.auto.T <- function(x) {
  glm1 <- glm(value ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_T , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(value ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}

###### CM_ext Temperate ####
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

## FD in Temperate ####
read.csv("data/spatial-survey-header-Tem.csv", sep= ";")%>%
  filter(plot%in%spatial_env_tem$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_T

glms.FD.T <- function(x) {
  glm1 <- glm(MPD ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw) + scale(disturbance), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_T , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(MPD ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(disturbance) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}
MPD.T(dark.dist.T, plot_x_sp_T)%>%
  rbind(MPD.T(WP.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(Tconstant.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(germtraits.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(seedmass.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(plantheight.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(leafarea.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(LDMC.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(SLA.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(Wplanttraits.dist.T, plot_x_sp_T))%>%
  rbind(MPD.T(alltraits.dist.T, plot_x_sp_T))%>%
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.T", "" , trait))%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, data_type, MPD, 
                elevation,FDD, GDD, Snw,  disturbance)%>%
  group_by(trait, data_type)%>%
  do(glms.FD.T (.))%>%
  write.csv("results/GLMs auto FD Tem.csv")
