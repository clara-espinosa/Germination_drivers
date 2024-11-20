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

# 5.2 Calculation of CM ####
sp_x_trait_M<- as.matrix(sp_x_trait_M)
plot_x_sp_M_PA<- as.matrix(plot_x_sp_M_PA)

CM_M<-functcomp(sp_x_trait_M, plot_x_sp_M_PA, CWM.type = "all") # CWM.type to consider binary or categorical traits
dim(plot_x_sp_M_PA)
# first let's check CM normality (all look normal distribution)
CM_M%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

# 5.3 Calculation of CWM ####
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



# 5.4 Dissimilarity matrices for MPD ####
# compute distance matrix separately for germination traits vs other plant traits
head(sp_x_trait_M)
sp_x_trait_M <- as.data.frame(sp_x_trait_M)
# individual traits distance matrices
dark.dist.M <- gowdis(sp_x_trait_M["odds_B_dark"])
WP.dist.M <- gowdis(sp_x_trait_M["odds_C_WP"])
Tconstant.dist.M <- gowdis(sp_x_trait_M["odds_D_constant"])

seedmass.dist.M <- gowdis(sp_x_trait_M["seed_mass"]) 
plantheight.dist.M <- gowdis(sp_x_trait_M["plant_height"]) 
leafarea.dist.M <- gowdis(sp_x_trait_M["leaf_area"]) 
LDMC.dist.M <- gowdis(sp_x_trait_M["LDMC"]) 
SLA.dist.M <- gowdis(sp_x_trait_M["SLA"]) 
# germination traits distance matrix
germtraits.dist.M <- gowdis(sp_x_trait_M[,6:8])

# plant traits weighted distance matrix by groups of traits 
# info from https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html
Wplanttraits.dist.M<-gawdis::gawdis(sp_x_trait_M[,1:5], w.type = "optimized", opti.maxiter = 200, 
                                    groups.weight=T, groups = c(1,2, 3, 3, 3))#

# 5.5. Calculation of MDP (main pairwise dissimilarity) each trait with melodic function ####
# Use of melodic function to compute both weighted and unweighted forms of MPD 
# activate function from melodic.R script in src folder
## 5.5. Calculation of MDP for individual traits
# try to do with a function
MPD.M <-function(distancematrix, communitydata){ # function to obtained omly interested indices from melodic function
  a <- deparse(substitute(distancematrix)) # get the name of the input distance matrix
  melodic(samp=communitydata, dis = distancematrix)-> x1 # plly melodic function to specific distance matrix
  as.data.frame(x1$abundance$mpd)%>%
    cbind(x1$presence$mpd)%>%
    mutate(trait = a)->x1
  colnames(x1) <- c("Abundance", "Presence", "trait")
  x1%>%
    cbind(plot_x_env_M2)->x1 # merge environmental data
  print(x1)
  
}

# would need to construct a loop feeding from a list of named matrices (work in progress)
# 1. list of matrices 
matrix.list.M <- list("Darkness"=dark.dist.M,"Water stress"= WP.dist.M)
for (i in matrix.list.M) {
  MPD.M(i,plot_x_sp_M )
}

# More manual but working!!
MPD.M(dark.dist.M, plot_x_sp_M)%>%
  rbind(MPD.M(WP.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(Tconstant.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(germtraits.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(seedmass.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(plantheight.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(leafarea.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(LDMC.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(SLA.dist.M, plot_x_sp_M))%>%
  rbind(MPD.M(Wplanttraits.dist.M, plot_x_sp_M))%>%
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.M", "" , trait))%>%
  mutate(trait= as.factor(trait))%>%
  mutate(trait= fct_recode(trait, "Germination traits"= "germtraits", "LDMC"= "LDMC", "SLA"="SLA", "Constant Temp"= "Tconstant",
                           "Water stress"="WP",  "Darkness"="dark","Leaf area"= "leafarea", 
                           "Plant height"= "plantheight", "Seed mass"= "seedmass", "Plant traits" = "Wplanttraits"))%>%
  mutate (trait = fct_relevel(trait, "Germination traits", "Darkness", "Water stress", "Constant Temp",
                              "Seed mass", "Plant height", "Leaf area", "LDMC", "SLA", "Plant traits"))-> FD.med
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
# 5.2.1 Calculation of CM  ####
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

# 5.4 Dissimilarity matrices for MPD ####
# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessary with our type of data
# compute distance matrix separately for germination traits vs other plant traits
head(sp_x_trait_T)
sp_x_trait_T <- as.data.frame(sp_x_trait_T)

# individual traits distance matrices
dark.dist.T <- gowdis(sp_x_trait_T["odds_B_dark"])
WP.dist.T <- gowdis(sp_x_trait_T["odds_C_WP"])
Tconstant.dist.T <- gowdis(sp_x_trait_T["odds_D_constant"])

seedmass.dist.T <- gowdis(sp_x_trait_T["seed_mass"]) 
plantheight.dist.T <- gowdis(sp_x_trait_T["plant_height"]) 
leafarea.dist.T <- gowdis(sp_x_trait_T["leaf_area"]) 
LDMC.dist.T <- gowdis(sp_x_trait_T["LDMC"]) 
SLA.dist.T <- gowdis(sp_x_trait_T["SLA"]) 

# germination traits distance matrix
germtraits.dist.T <- gowdis(sp_x_trait_T[,6:8])

# plant traits weighted distance matrix by groups of traits 
# info from https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html
Wplanttraits.dist.T<-gawdis(sp_x_trait_T[,1:5], w.type = "optimized", opti.maxiter = 200, 
                            groups.weight=T, groups = c(1,2, 3, 3, 3))

# all traits distance matrix 
alltraits.dist.T <- gawdis::gawdis(as.data.frame(sp_x_trait_T[,1:8]), w.type = "optimized", opti.maxiter = 200, 
                                   groups.weight=T, groups = c(1,2, 3, 3, 3, 4, 5, 6))


# 5.5. Calculation of MDP (main pairwise dissimilarity) each trait with melodic function ####
# Use of melodic function to compute both weighted and no weighted forms of MPD and Rao
# activate function from melodic.R script in src folder
## 5.5. Calculation of MDP for individual traits
# try to do with a function

MPD.T <-function(distancematrix, communitydata){ # function to obtained omly interested indices from melodic function
  a <- deparse(substitute(distancematrix)) # get the name of the input distance matrix
  melodic(samp=communitydata, dis = distancematrix)-> x1 # plly melodic function to specific distance matrix
  as.data.frame(x1$abundance$mpd)%>%
    cbind(x1$presence$mpd)%>%
    mutate(trait = a)->x1
  colnames(x1) <- c("Abundance", "Presence", "trait")
  x1%>%
    cbind(plot_x_env_T2)->x1 # merge environmental data
  print(x1)
  
}

# would need to construct a loop feeding from a list of named matrices (work in progress)
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
  mutate(trait= as.factor(trait))%>%
  mutate(trait= fct_recode(trait, "Germination traits"= "germtraits", "LDMC"= "LDMC", "SLA"="SLA", "Constant Temp"= "Tconstant",
                           "Water stress"="WP",  "Darkness"="dark","Leaf area"= "leafarea", 
                           "Plant height"= "plantheight", "Seed mass"= "seedmass", "Plant traits" = "Wplanttraits"))%>%
  mutate (trait = fct_relevel(trait, "Germination traits", "Darkness", "Water stress", "Constant Temp",
                              "Seed mass", "Plant height", "Leaf area", "LDMC", "SLA", "Plant traits"))->FD.tem
