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
# 5.4 Dissimilarity matrices for MPD and RaoQ ####
# compute distance matrix separately for germination traits vs other plant traits
head(sp_x_trait_M)
sp_x_trait_M <- as.data.frame(sp_x_trait_M)

# germination traits distance matrix
germ.traits.dist.M <- gowdis(sp_x_trait_M[,6:8])

# plant traits NO weighted distance matrix
plant.traits.dist.M <- gowdis(sp_x_trait_M[,1:5]) 

# plant traits weighted distance matrix by groups of traits 
# info from https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html
plant.traits.weighted.dist.M<-gawdis(sp_x_trait_M[,1:5], w.type = "optimized", opti.maxiter = 200, 
                                     groups.weight=T, groups = c(1,2, 3, 3, 3))#

# 5.5. Calculation of MDP (main pairwise dissimilarity) and RAo with melodic function ####
# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
## 5.5.1 Germ traits####
plot_x_sp_M <- as.data.frame(plot_x_sp_M)
germ.melodic.M <- melodic(plot_x_sp_M, germ.traits.dist.M)
str(germ.melodic.M)

# germ melodic object data handling
as.data.frame(germ.melodic.M$abundance$rao) %>%
  cbind (as.data.frame(germ.melodic.M$abundance$mpd))%>%
  #cbind (as.data.frame(germ.melodic.M$presence$rao))%>%
  #cbind (as.data.frame(germ.melodic.M$presence$mpd))%>%
  cbind (plot_x_env_M2)%>%
  mutate(Germ.Rao.ab = germ.melodic.M$abundance$rao)%>%
  mutate(Germ.Mpd.ab = germ.melodic.M$abundance$mpd)%>%
  #mutate(Rao.pa = germ.melodic.M$abundance$mpd)%>%
  #mutate(Mpd.pa = germ.melodic.M$abundance$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Germ.Rao.ab,Germ.Mpd.ab)%>% #, Rao.pa, Mpd.pa
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot,Germ.Rao.ab, Germ.Mpd.ab,  elevation, FDD, GDD, Snw)-> FD.germ.M

## 5.5.2 Plant traits Weighted  ####
weighted.plant.melodic.M <- melodic(plot_x_sp_M, plant.traits.weighted.dist.M)
str(weighted.plant.melodic.M)

# plant melodic object data handling
as.data.frame(weighted.plant.melodic.M$abundance$rao) %>%
  cbind (as.data.frame(weighted.plant.melodic.M$abundance$mpd))%>%
  cbind (plot_x_env_M2)%>%
  mutate(W.plant.Rao.ab = weighted.plant.melodic.M$abundance$rao)%>%
  mutate(W.plant.Mpd.ab = weighted.plant.melodic.M$abundance$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,W.plant.Rao.ab,W.plant.Mpd.ab)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot,W.plant.Rao.ab,W.plant.Mpd.ab, elevation, FDD, GDD, Snw)-> FD.w.plant.M

## Join germ and plant mdp and rao calculations
FD.germ.M %>%
  merge(FD.w.plant.M, by = c("site", "plot", "elevation", "FDD","GDD", "Snw"))-> FD.M

# first let's check MPD and Rao normality
x11()
FD.M%>%
  gather(Func.div, value, Germ.Rao.ab:W.plant.Mpd.ab)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = Func.div), color = "black") + 
  facet_wrap(~Func.div, ncol =2, scales = "free") +
  theme_bw(base_size = 12) # somewhat normal distributed !?!?!

# data handling for visualization Rao and MPD (weigthed by abundance) vs microclimatic gradients.
ann.sig.FD.M <- data.frame (read.csv("results/sig_FD_M.csv", sep = ";"))

FD.M%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Rao.ab:W.plant.Mpd.ab)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(Func.div~micro_variable, scales = "free")+
  geom_text(data=ann.sig.FD.M, label =ann.sig.FD.M$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Mpd and Rao in Mediterranean system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_blank(),
        legend.position = "bottom")-> FD_micro_M;FD_micro_M

ggsave(FD_micro_M, file = "FD vs micro Med.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
lms.fd <- function(x) {
  glm(value ~ elevation+ FDD + GDD + Snw, data = x) -> m1
  broom::tidy(m1)
}
FD.M%>%
  gather (Func.div, value, Germ.Rao.ab:W.plant.Mpd.ab)%>%
  group_by(Func.div)%>%
  do(lms.fd(.)) %>%
  write.csv("results/FD vs micro MED.csv")

## 5.5.3 Plant traits NO weighted (not use?) ##### 
plot_x_sp_M <- as.matrix(plot_x_sp_M)
plant.melodic.M <- melodic(plot_x_sp_M, plant.traits.dist.M)
str(plant.melodic.M)

as.data.frame(plant.melodic.M$abundance$rao) %>%
  cbind (as.data.frame(plant.melodic.M$abundance$mpd))%>%
  cbind (plot_x_env_M2)%>%
  mutate(Rao.ab = plant.melodic.M$abundance$rao)%>%
  mutate(Mpd.ab = plant.melodic.M$abundance$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Rao.ab,Mpd.ab)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot,Rao.ab,Mpd.ab, elevation, FDD, GDD, Snw)%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Rao.ab:Mpd.ab)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(Func.div~micro_variable, scales = "free")+
  #geom_text(data=ann.sig.CWM.T, label =ann.sig.CWM.T$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Mpd and Rao in Mediterranean system for unweighted plant traits")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")


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
# 5.4 Dissimilarity matrices for MPD and RaoQ ####
# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessary with our type of data
# compute distance matrix separately for germination traits vs other plant traits
head(sp_x_trait_T)
sp_x_trait_T <- as.data.frame(sp_x_trait_T)

# germination traits distance matrix
germ.traits.dist.T <- gowdis(sp_x_trait_T[,6:8])

# plant traits NO weighted distance matrix
plant.traits.dist.T <- gowdis(sp_x_trait_T[,1:5]) 

# plant traits weighted distance matrix by groups of traits 
# info from https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html
plant.traits.weighted.dist.T<-gawdis(sp_x_trait_T[,1:5], w.type = "optimized", opti.maxiter = 200, 
                                     groups.weight=T, groups = c(1,2, 3, 3, 3))

# 5.5. Calculation of MDP (main pairwise dissimilarity) and RAo with melodic function ####
# Use of melodic function to compute both weighted and no weighted forms of MPD and Rao
# activate function from melodic.R script in src folder
## 5.5.1 Germ traits####
plot_x_sp_T <- as.data.frame(plot_x_sp_T)
germ.melodic.T <- melodic(plot_x_sp_T, germ.traits.dist.T)
str(germ.melodic.T)

# germ melodic object data handling
as.data.frame(germ.melodic.T$abundance$rao) %>%
  cbind (as.data.frame(germ.melodic.T$abundance$mpd))%>%
  #cbind (as.data.frame(germ.melodic.T$presence$rao))%>%
  #cbind (as.data.frame(germ.melodic.T$presence$mpd))%>%
  cbind (plot_x_env_T2)%>%
  mutate(Germ.Rao.ab = germ.melodic.T$abundance$rao)%>%
  mutate(Germ.Mpd.ab = germ.melodic.T$abundance$mpd)%>%
  #mutate(Rao.pa = germ.melodic.T$abundance$mpd)%>%
  #mutate(Mpd.pa = germ.melodic.T$abundance$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Germ.Rao.ab,Germ.Mpd.ab)%>% #, Rao.pa, Mpd.pa
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot,Germ.Rao.ab, Germ.Mpd.ab,  elevation, FDD, GDD, Snw)-> FD.germ.T

## 5.5.2 Plant traits Weighted  ####
weighted.plant.melodic.T <- melodic(plot_x_sp_T, plant.traits.weighted.dist.T)
str(weighted.plant.melodic.T)

# plant melodic object data handling
as.data.frame(weighted.plant.melodic.T$abundance$rao) %>%
  cbind (as.data.frame(weighted.plant.melodic.T$abundance$mpd))%>%
  cbind (plot_x_env_T2)%>%
  mutate(W.plant.Rao.ab = weighted.plant.melodic.T$abundance$rao)%>%
  mutate(W.plant.Mpd.ab = weighted.plant.melodic.T$abundance$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,W.plant.Rao.ab,W.plant.Mpd.ab)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot,W.plant.Rao.ab,W.plant.Mpd.ab, elevation, FDD, GDD, Snw)-> FD.w.plant.T

## Join germ and plant mdp and rao calculations
FD.germ.T %>%
  merge(FD.w.plant.T, by = c("site", "plot", "elevation", "FDD","GDD", "Snw"))-> FD.T

# first let's check MPD and Rao normality
x11()
FD.T%>%
  gather(Func.div, value, Germ.Rao.ab:W.plant.Mpd.ab)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = Func.div), color = "black") + 
  facet_wrap(~Func.div, ncol =2, scales = "free") +
  theme_bw(base_size = 12) #  normally distributed 

# data handling for visualization Rao and MPD (weigthed by abundance) vs microclimatic gradients.
ann.sig.FD.T <- data.frame (read.csv("results/sig_FD_T.csv", sep = ";"))

FD.T%>%
  mutate(site = as.factor(site))%>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Rao.ab:W.plant.Mpd.ab)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(Func.div~micro_variable, scales = "free")+
  geom_text(data=ann.sig.FD.T, label =ann.sig.FD.T$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Mpd and Rao in Temperate system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_blank(),
        legend.position = "bottom")-> FD_micro_T;FD_micro_T

ggsave(FD_micro_T, file = "FD vs micro Tem.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
lms.fd <- function(x) {
  glm(value ~ elevation+ FDD + GDD + Snw, data = x) -> m1
  broom::tidy(m1)
}
FD.T%>%
  gather (Func.div, value, Germ.Rao.ab:W.plant.Mpd.ab)%>%
  group_by(Func.div)%>%
  do(lms.fd(.)) %>%
  write.csv("results/FD vs micro TEM.csv")

## 5.5.3 Plant traits NO weighted (not use?) ####
plot_x_sp_T <- as.data.frame(plot_x_sp_T)
plant.melodic.T <- melodic(plot_x_sp_T, plant.traits.dist.T)
str(plant.melodic.T)

as.data.frame(plant.melodic.T$abundance$rao) %>%
  cbind (as.data.frame(plant.melodic.T$abundance$mpd))%>%
  cbind (plot_x_env_T2)%>%
  mutate(Rao.ab = plant.melodic.T$abundance$rao)%>%
  mutate(Mpd.ab = plant.melodic.T$abundance$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Rao.ab,Mpd.ab)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot,Rao.ab,Mpd.ab, elevation, FDD, GDD, Snw)%>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Rao.ab:Mpd.ab)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(Func.div~micro_variable, scales = "free")+
  #geom_text(data=ann.sig.CWM.T, label =ann.sig.CWM.T$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Mpd and Rao in Temperate system for unweighted plant traits")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")

