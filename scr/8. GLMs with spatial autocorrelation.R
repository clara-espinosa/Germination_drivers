library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis); library(glmm); library (lme4); library(spdep) # para crear la autocovariate
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
# FD apart due to 2 extra additional traits  (germ vs plant) as explanatory variable
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
  group_by(trait,community_metric)%>%
  do(glms.auto.M(.))%>%
  write.csv("results/GLMs auto Med.csv")
  

glms.auto.M <- function(x) {
  glm1 <- glm(value ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_M , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(value ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
  broom::tidy(glm1.auto)
}


## FD in Mediterranean ####
read.csv("data/spatial-survey-header-Med.csv")%>%
  filter(plot%in%spatial_env_med$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_M

glms.FD.M <- function(x) {
  glm1 <- glm(MPD ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_M , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(MPD ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
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
  #rbind(MPD.M(alltraits.dist.M, plot_x_sp_M))%>%
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.M", "" , trait))%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, data_type, MPD, 
                elevation,FDD, GDD, Snw)%>%
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

## FD in Temperate ####
read.csv("data/spatial-survey-header-Tem.csv", sep= ";")%>%
  filter(plot%in%spatial_env_tem$plot)%>% 
  dplyr::select(longitude, latitude)%>% 
  as.matrix ()->longlat_T

glms.FD.T <- function(x) {
  glm1 <- glm(MPD ~ scale(elevation) +scale(FDD) + scale(GDD) +  scale(Snw), data = x) 
  autocov <- autocov_dist(glm1$residuals, longlat_T , nbs = 5, type= "inverse", style = "B", zero.policy = T, longlat = TRUE)
  x$auto <- autocov
  glm1.auto<- glm(MPD ~ scale(elevation) + scale(FDD) + scale(GDD) +  scale(Snw) + scale(auto), data = x)
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
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.T", "" , trait))%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, latitude, longitude,
                trait, data_type, MPD, 
                elevation,FDD, GDD, Snw)%>%
  group_by(trait, data_type)%>%
  do(glms.FD.T (.))%>%
  write.csv("results/GLMs auto FD Tem.csv")


########## GLMs SIZE EFFECTS VISUALIZATION ######
### Community metrics Mediterranean ####
effect_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                  "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                  "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                  "LDMC" = "LDMC", "SLA" = "SLA", "CM_ext"= "Community Means Extended",
                  "CM" = "Communty Means", "CWM"="Community Weighted Means")
x11()
read.csv("results/GLMs auto Med.csv", sep =",")%>%
  mutate(data_type = ifelse(community_metric== "CM", "Presence", "Abundance"))%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  mutate(community_metric = factor(community_metric))%>%
  mutate(community_metric = fct_relevel(community_metric, "CM", "CWM"))%>%
  filter(!term == "(Intercept)")%>%
  filter(!term == "scale(auto)")%>%
  mutate(term= factor(term))%>%
  mutate(term= fct_recode(term, "Elevation" = "scale(elevation)", "FDD"="scale(FDD)",
                      "GDD"="scale(GDD)", "Snow"="scale(Snw)"))%>%
  mutate(CI = 1.96*std.error)%>% # confidence interval multiply 1.96 per std error
  mutate(CImin = estimate-CI, 
         CImax= estimate+CI)%>%
  mutate(color = case_when(CImin>0 & CImax>0 ~ "deepskyblue",
                           CImin<0 & CImax<0~ "deepskyblue",
                           TRUE ~"grey"))%>%
  ggplot(aes(x= term, y =estimate, ymin = CImin, ymax = CImax))+
  geom_errorbar (aes(color=color),width = 0, linewidth =1.2) + #, color="black" 
  geom_point(aes(color=color), size = 3) +#
  scale_color_manual (values = c("deepskyblue"="deepskyblue","deepskyblue"="deepskyblue", "grey"= "grey" ))+
  facet_grid (community_metric~trait,labeller = as_labeller(effect_names),  scales = "free_x") + #ncol = 8, nrow =2,
  geom_hline(yintercept = 0, linetype = "dashed", linewidth =1, color = "red") +
  coord_flip() +
  labs(y = "Parameter estimate", title = "Mediterranean community") + # Temperate   Mediterranean
  theme_classic (base_size = 12)+#ggthemes::theme_tufte(base_size = 14) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size = 20),
        strip.text = element_text( size = 16), #face = "bold",
        strip.text.y = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank (),
        panel.background = element_rect(color = "black", fill = NULL),
        plot.tag.position = c(0.015,1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", angle = 30),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text (size=14))

### Community metrics Temperate ####
effect_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                  "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                  "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                  "LDMC" = "LDMC", "SLA" = "SLA", "CM_ext"= "Community Means Extended",
                  "CM" = "Communty Means", "CWM"="Community Weighted Means")
x11()
read.csv("results/GLMs auto Tem.csv", sep =";")%>%
  mutate(data_type = ifelse(community_metric== "CM", "Presence", "Abundance"))%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  mutate(community_metric = factor(community_metric))%>%
  mutate(community_metric = fct_relevel(community_metric, "CM", "CWM"))%>%
  filter(!term == "(Intercept)")%>%
  filter(!term == "scale(auto)")%>%
  mutate(term= factor(term))%>%
  mutate(term= fct_recode(term, "Elevation" = "scale(elevation)", "FDD"="scale(FDD)",
                          "GDD"="scale(GDD)", "Snow"="scale(Snw)"))%>%
  mutate(CI = 1.96*std.error)%>% # confidence interval multiply 1.96 per std error
  mutate(CImin = estimate-CI, 
         CImax= estimate+CI)%>%
  mutate(color = case_when(CImin>0 & CImax>0 ~ "deepskyblue",
                           CImin<0 & CImax<0~ "deepskyblue",
                           TRUE ~"grey"))%>%
  ggplot(aes(x= term, y =estimate, ymin = CImin, ymax = CImax))+
  geom_errorbar (aes(color=color),width = 0, linewidth =1.2) + #, color="black" 
  geom_point(aes(color=color), size = 3) +#
  scale_color_manual (values = c("deepskyblue"="deepskyblue","deepskyblue"="deepskyblue", "grey"= "grey" ))+
  facet_grid (community_metric~trait,labeller = as_labeller(effect_names),  scales = "free_x") + #ncol = 8, nrow =2,
  geom_hline(yintercept = 0, linetype = "dashed", linewidth =1, color = "red") +
  coord_flip() +
  labs(y = "Parameter estimate", title = "Temperate community") + # Temperate   Mediterranean
  theme_classic (base_size = 12)+#ggthemes::theme_tufte(base_size = 14) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size = 20),
        strip.text = element_text( size = 16), #face = "bold",
        strip.text.y = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank (),
        panel.background = element_rect(color = "black", fill = NULL),
        plot.tag.position = c(0.015,1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", angle = 30),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text (size=14))


### FD Presence - abundance MEDITERRANEAN ####
x11()
read.csv("results/GLMs auto FD Med.csv", sep=",") %>% # table with glm model results from script 7
  convert_as_factor(trait, data_type,term) %>%
  filter(!term == "(Intercept)")%>%
  filter(!term == "scale(auto)")%>%
  mutate(trait= fct_recode(trait, "LDMC"= "LDMC", "SLA"="SLA", "Constant Temp"= "Tconstant",
                           "Water stress"="WP", "Adult plant"="Wplanttraits", "Darkness"="dark",
                           "Germ traits" = "germtraits", "Leaf area"= "leafarea", 
                           "Plant height"= "plantheight", "Seed mass"= "seedmass"))%>%
  mutate (trait = fct_relevel(trait, "Germ traits", "Darkness", "Water stress", "Constant Temp",
                              "Adult plant", "Seed mass", "Plant height", "Leaf area", "LDMC", "SLA"))%>%
  mutate(term= fct_recode(term, "Elevation" = "scale(elevation)", "FDD"="scale(FDD)",
                          "GDD"="scale(GDD)", "Snow"="scale(Snw)"))%>%
  mutate(CI = 1.96*std.error)%>% # confidence interval multiply 1.96 per std error
  mutate(CImin = estimate-CI, 
         CImax= estimate+CI)%>%
  mutate(color = case_when(CImin>0 & CImax>0 ~ "deepskyblue",
                           CImin<0 & CImax<0~ "deepskyblue",
                           TRUE ~"grey"))%>%
  ggplot(aes(x= term, y =estimate, ymin = CImin, ymax = CImax))+
  geom_errorbar (aes(color=color),width = 0, linewidth =1.2) + #, color="black" 
  geom_point(aes(color=color), size = 3) +#
  scale_color_manual (values = c("deepskyblue"="deepskyblue","deepskyblue"="deepskyblue", "grey"= "grey" ))+
  facet_grid (data_type~trait,  scales = "free_x") + #ncol = 8, nrow =2,
  geom_hline(yintercept = 0, linetype = "dashed", linewidth =1, color = "red") +
  coord_flip() +
  labs(y = "Parameter estimate", title = "Mediterranean community") + # Temperate   Mediterranean
  theme_classic (base_size = 12)+#ggthemes::theme_tufte(base_size = 14) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size = 20),
        strip.text = element_text( size = 16), #face = "bold",
        strip.text.y = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank (),
        panel.background = element_rect(color = "black", fill = NULL),
        plot.tag.position = c(0.015,1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", angle = 30),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text (size=14))

### FD Presence - abundance Temperate ####
x11()
read.csv("results/GLMs auto FD Tem.csv", sep=",") %>% # table with glm model results from script 7
  convert_as_factor(trait, data_type,term) %>%
  filter(!term == "(Intercept)")%>%
  filter(!term == "scale(auto)")%>%
  filter(!trait == "alltraits")%>%
  mutate(trait= fct_recode(trait, "LDMC"= "LDMC", "SLA"="SLA", "Constant Temp"= "Tconstant",
                           "Water stress"="WP", "Adult plant"="Wplanttraits", "Darkness"="dark",
                           "Germ traits" = "germtraits", "Leaf area"= "leafarea", 
                           "Plant height"= "plantheight", "Seed mass"= "seedmass"))%>%
  mutate (trait = fct_relevel(trait, "Germ traits", "Darkness", "Water stress", "Constant Temp",
                              "Adult plant", "Seed mass", "Plant height", "Leaf area", "LDMC", "SLA"))%>%
  mutate(term= fct_recode(term, "Elevation" = "scale(elevation)", "FDD"="scale(FDD)",
                          "GDD"="scale(GDD)", "Snow"="scale(Snw)"))%>%
  mutate(CI = 1.96*std.error)%>% # confidence interval multiply 1.96 per std error
  mutate(CImin = estimate-CI, 
         CImax= estimate+CI)%>%
  mutate(color = case_when(CImin>0 & CImax>0 ~ "deepskyblue",
                           CImin<0 & CImax<0~ "deepskyblue",
                           TRUE ~"grey"))%>%
  ggplot(aes(x= term, y =estimate, ymin = CImin, ymax = CImax))+
  geom_errorbar (aes(color=color),width = 0, linewidth =1.2) + #, color="black" 
  geom_point(aes(color=color), size = 3) +#
  scale_color_manual (values = c("deepskyblue"="deepskyblue","deepskyblue"="deepskyblue", "grey"= "grey" ))+
  facet_grid (data_type~trait,  scales = "free_x") + #ncol = 8, nrow =2,
  geom_hline(yintercept = 0, linetype = "dashed", linewidth =1, color = "red") +
  coord_flip() +
  labs(y = "Parameter estimate", title = "Temperate community") + # Temperate   Mediterranean
  theme_classic (base_size = 12)+#ggthemes::theme_tufte(base_size = 14) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size = 20),
        strip.text = element_text( size = 16), #face = "bold",
        strip.text.y = element_text(size = 14),
        legend.position = "none",
        strip.background = element_blank (),
        panel.background = element_rect(color = "black", fill = NULL),
        plot.tag.position = c(0.015,1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black", angle = 30),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text (size=14))
