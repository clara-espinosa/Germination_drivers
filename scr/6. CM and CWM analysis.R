library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis)
###################################  MEDITERRANEAN ###############################################################
# Following chapter 5 R TRait-Based ecology by Lars####
# 5.1 Data matrices ####
# in sp x plot 
# in all files names of species and plots are not part of the matrix, used as row and columns names
# mediterranean #
sp_x_plot_M # species in rows and plots in columns, abundance already as relative cover
plot_x_sp_M  # plot in row and species in columns, abundance in relative cover data
plot_x_sp_M_PA # plot in rw and species in columns, presence-abscence data
sp_x_trait_M # species in rows and traits in columns all traits already transformed
# remove species which had 0 germination across all treatments in germination experiment Solidago virgaurea
plot_x_env_M # plot in rows and env data in columns

# extended dataset for CM
sp_x_trait_M
plot_x_env_M_ext
plot_x_sp_M_PA_ext

## 5.1.1 Check normality of the quantitative traits ####
### all traits are right skewed, try a log transformation look better but still not normal
sp_x_trait_M%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# odds do not look normal distributed!!!

# 5.2.1 Calculation of CM restricted dataset ####
sp_x_trait_M<- as.matrix(sp_x_trait_M)
plot_x_sp_M_PA<- as.matrix(plot_x_sp_M_PA)

CM_M<-functcomp(sp_x_trait_M, plot_x_sp_M_PA, CWM.type = "all") # CWM.type to consider binary or categorical traits

# first let's check CM normality (all look normal distribution)
CM_M%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# all traits look quite normally distributed!!
plot_x_env_M%>%
  rownames_to_column(var = "plot")->plot_x_env_M2 

# loop for visualization plots
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

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model 
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cm <- function(x) {
  glm(CM ~ elevation +FDD + GDD +  Snw, data = x) -> m1
  broom::tidy(m1)
}

CM_M%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_M2, by= "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cm(.)) %>%
  write.csv("results/CM vs micro MED.csv")

# double facetting with significances
ann.sig.CM.M <- data.frame (read.csv("results/sig_CM_M.csv", sep = ";"))
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                  "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")
str(CM_microclima_M)
x11()
CM_microclima_M%>%
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

rda_med_all<-rda(CM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
plot(rda_med_all, type = "n", scaling = "sites")
text(rda_med_all, dis = "cn", scaling = "sites")
text(rda_med_all, dis = "sp", scaling = "sites", col = "red")

rda_med_0<-rda(CM_M~1, data= plot_x_env_M)
rda_med_all<-rda(CM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
ordistep(rda_med_0, scope=formula(rda_med_all), direction = "forward")

RsquareAdj (rda_med_all)$adj.r.squared # 0.06
# final model only with elevation

# 5.2.2 Calculation of CM extended dataset ####
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

## 5.3 Calculation of CWM ####
##dataframes in matrix format
sp_x_trait_M <- as.matrix(sp_x_trait_M)
plot_x_sp_M <- as.matrix(plot_x_sp_M)
CWM_M<-functcomp(sp_x_trait_M, plot_x_sp_M, CWM.type = "all") # CWM.type to consider binary or categorical traits

# first let's check CWM normality
x11()
CWM_M%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# odds do not look normally distributed
#PROBLEM?!?!? 

# Let’s see if CWM changes along the climatic gradient
# explore visually a bit the results before running any analysis,
plot_x_env_M%>%
  rownames_to_column(var = "plot")->plot_x_env_M2 

# loop for visualization plots
CWM_M%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_M2)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "plant_height", "seed_mass",
                             "leaf_area", "LDMC", "SLA",
                             "odds_B_dark","odds_C_WP","odds_D_constant"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CWM_microclima_M

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cwm <- function(x) {
  glm(CWM ~ elevation+ FDD + GDD + Snw, data = x) -> m1
  broom::tidy(m1)
}

CWM_M%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_M2, by= "plot")%>%
  gather(trait, CWM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cwm(.)) %>%
  write.csv("results/CWM vs micro MED.csv")

# see summary in results
# double facetting with significances
ann.sig.CWM.M <- data.frame (read.csv("results/sig_CWM_M.csv", sep = ";"))
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                 "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")
CWM_microclima_M%>%
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
  geom_text(data=ann.sig.CWM.M, label =ann.sig.CWM.M$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Community Weighted Means in Mediterranean system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> plot_CWM_M;plot_CWM_M
ggsave(plot_CWM_M, file = "CWM Med trait vs micro.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
#width = 320, height = 260, units = "mm", 

# we can also use the multivariate analyses for the representation, for example RDA

rda_med_all<-rda(CWM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
plot(rda_med_all, type = "n", scaling = "sites")
text(rda_med_all, dis = "cn", scaling = "sites")
text(rda_med_all, dis = "sp", scaling = "sites", col = "red")

rda_med_0<-rda(CWM_M~1, data= plot_x_env_M)
rda_med_all<-rda(CWM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
ordistep(rda_med_0, scope=formula(rda_med_all), direction = "forward")

RsquareAdj (rda_med_all)$adj.r.squared # 0.38
# final model with elevation, GDD and Snow significant and thus included!!

###################################  TEMPERATE  ##################################################################
# Following chapter 5 R TRait-Based ecology by Lars####
# 5.1 Data matrices ####
# base matrices from species level analysis but with reformatting
# in sp x plot (species in rows and plots in columns)
# in all files names of species and plots are not part of the matrix, used as row and columns names

sp_x_plot_T# species in rows and plots in columns, abundance already as relative cover
plot_x_sp_T# plot in row and species in columns, abundance in relative cover data
plot_x_sp_T_PA # plot in row and species in columns, presence-absence data
sp_x_trait_T# species in rows and traits in columns remove species which had 0 germination across al treatments in germination experiment
plot_x_env_T# plot in rows and env data in columns


# extended dataset for CM
sp_x_trait_T
plot_x_env_T_ext
plot_x_sp_T_PA_ext

## 5.1.1 Check normality of the quantitative traits ####
### all traits are right skewed, ty a log transformation look better but still not normal
sp_x_trait_T%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

# odds do not look normally distributed
# 5.2.1 Calculation of CM  restricted dataset ####
sp_x_trait_T<- as.matrix(sp_x_trait_T)
plot_x_sp_T_PA<- as.matrix(plot_x_sp_T_PA)

CM_T<-functcomp(sp_x_trait_T, plot_x_sp_T_PA, CWM.type = "all") # CWM.type to consider binary or categorical traits

# first let's check CWM normality (all look normal distribution)
CM_T%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# all look pretty normal distributed

plot_x_env_T%>%
  rownames_to_column(var = "plot")->plot_x_env_T2 

CM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra")) %>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait,"odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass", "plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CM_microclima_T

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cm <- function(x) {
  glm(CM ~ elevation+ FDD + GDD+ Snw, data = x) -> m1
  broom::tidy(m1)
}

CM_T%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_T2, by= "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cm(.)) %>%
  write.csv("results/CM vs micro TEM.csv")
# double facetting with significances
ann.sig.CM.T <- data.frame (read.csv("results/sig_CM_T.csv", sep = ";"))
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                 "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")

CM_microclima_T%>%
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

rda_tem_all<-rda(CM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
plot(rda_tem_all, type = "n", scaling = "sites")
text(rda_tem_all, dis = "cn", scaling = "sites")
text(rda_tem_all, dis = "sp", scaling = "sites", col = "red")

rda_tem_0<-rda(CM_T~1, data= plot_x_env_T)
rda_tem_all<-rda(CM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
ordistep(rda_tem_0, scope=formula(rda_tem_all), direction = "forward")

RsquareAdj (rda_tem_all)$adj.r.squared # 0.26
# final model with GDD and elevation significant and thus included!!
# 5.2.2 Calculation of CM  extended dataset ####
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
# 5.3 Calculation of CWM ####
sp_x_trait_T<- as.matrix(sp_x_trait_T)
plot_x_sp_T <- as.matrix(plot_x_sp_T)
CWM_T<-functcomp(sp_x_trait_T, plot_x_sp_T, CWM.type = "all") # CWM.type to consider binary or categorical traits

# first let's check CWM normality (all look normal distribution)
x11()
CWM_T%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# generally quite normal distributed

# Let’s see if CWM changes along the climatic gradient
# explore visually a bit the results before running any analysis,
plot_x_env_T%>%
  rownames_to_column(var = "plot")->plot_x_env_T2 

CWM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra")) %>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait,"odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass", "plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CWM_microclima_T
# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model

lms.cwm <- function(x) {
  glm(CWM ~ elevation+ FDD + GDD + Snw, data = x) -> m1
  broom::tidy(m1)
}

CWM_T%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_T2, by= "plot")%>%
  gather(trait, CWM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cwm(.)) %>%
  write.csv("results/CWM vs micro TEM.csv")
# see summary in results

# double facetting with significances
ann.sig.CWM.T <- data.frame (read.csv("results/sig_CWM_T.csv", sep = ";"))
trait_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                 "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                 "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                 "LDMC" = "LDMC", "SLA" = "SLA", "elevation"= "Elevation", 
                 "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")

CWM_microclima_T%>%
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
  geom_text(data=ann.sig.CWM.T, label =ann.sig.CWM.T$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Community Weighted Means in Temperate system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> plot_CWM_T;plot_CWM_T
ggsave(plot_CWM_T, file = "CWM Tem trait vs micro.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 

# we can also use the multivariate analyses for the representation, for example RDA

rda_tem_all<-rda(CWM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
plot(rda_tem_all, type = "n", scaling = "sites")
text(rda_tem_all, dis = "cn", scaling = "sites")
text(rda_tem_all, dis = "sp", scaling = "sites", col = "red")

rda_tem_0<-rda(CWM_T~1, data= plot_x_env_T)
rda_tem_all<-rda(CWM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
ordistep(rda_tem_0, scope=formula(rda_tem_all), direction = "forward")

RsquareAdj (rda_tem_all)$adj.r.squared # 0.01
# with GDD is the best model!



# Effect size plots for CM and CWM ####
effect_names <- c("odds_B_dark" = "Darkness","odds_C_WP" = "Water stress",
                  "odds_D_constant" = "Constant Temp","seed_mass" = "Seed mass",
                  "plant_height" = "Plant height", "leaf_area" = "Leaf area",
                  "LDMC" = "LDMC", "SLA" = "SLA", "CM" = "Communty Means", "CWM"="Community Weighted Means")
dat_text <- data.frame(label = "positive", "negative")

x11()
read.csv("results/lm community estimates graphs.csv", sep=";") %>% # table with lm model results from script 6 CM 
  convert_as_factor(community, analysis, trait, term) %>%
  mutate(trait = fct_relevel(trait, "odds_B_dark","odds_C_WP","odds_D_constant",
                             "seed_mass","plant_height", 
                             "leaf_area", "LDMC", "SLA"))%>%
  mutate(term=fct_recode(term, "Elevation"="elevation", 
                         "FDD" = "FDD", "GDD"="GDD", "Snow"="Snw"))%>%
  filter(community == "Mediterranean")%>%
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
  facet_grid (analysis~trait,labeller = as_labeller(effect_names),  scales = "free_x") + #ncol = 8, nrow =2,
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
        axis.title.x = element_text (size=14))-> Comsig_M;Comsig_M #Comsig_T;Comsig_T

ggsave(Comsig_M, file = "Mediterranean community metrics significances.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 