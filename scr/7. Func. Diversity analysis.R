library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis):library(purrr);library(stringr):library(tibble)
###################################  MEDITERRANEAN ###############################################################
# Following chapter 5 R Trait-Based ecology by Lars####
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
as.data.frame(sp_x_trait_M)%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)
# odds do not look normal distributed!!!
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
rm(all.traits.dist.M)
# all traits distance matrix 
alltraits.dist.M <- gawdis::gawdis(as.data.frame(sp_x_trait_M[,1:8]), w.type = "optimized", opti.maxiter = 200, 
                                    groups.weight=T, groups = c(1,2, 3, 3, 3, 4, 5, 6))

# check distribution of distances germ vs plant (not necessary?)
as.matrix(germ.traits.dist.M)%>%
  as.data.frame(germ.traits.dist.M)%>%
  gather(sp, germ.dist, 1:18)%>%
  cbind((as.matrix(plant.traits.weighted.dist.M)%>%
          as.data.frame(plant.traits.weighted.dist.M)%>%
          gather(sp, plant.dist, 1:18)))%>%
  dplyr::select(1,2,4)%>%
  gather(trait, dist, germ.dist:plant.dist)%>%
  filter(dist>0)%>%
  ggplot(aes(dist, fill=trait))+
  #geom_density(alpha= 0.5)+
  geom_histogram(color = "black")+
  scale_fill_manual(values = c("dodgerblue4","limegreen"))

# 5.5. Calculation of MDP (main pairwise dissimilarity) with melodic function ####
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
  rbind(MPD.M(alltraits.dist.M, plot_x_sp_M))%>%
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.M", "" , trait))%>%
  group_by(trait, data_type)%>%
  get_summary_stats(MPD)%>%
  mutate(trait= as.factor(trait))%>%
  mutate(trait = fct_relevel(trait, "alltraits", "Wplanttraits", "germtraits",
                             "dark", "WP","Tconstant", 
                             "seedmass", "plantheight", "leafarea", "LDMC", "SLA"))%>%
  ggplot()+
  geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_point(aes(x= trait, y= mean, fill = trait), color = "black", shape= 21, size = 5)+
  #scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  scale_y_continuous(limits = c(0.05, 0.7))+
  coord_flip()+
  facet_grid(~data_type)+
  labs(title = "Mediterranean community FD differences", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        #axis.title.x = element_blank(),
        legend.position = "none")-> Full.FD.M;Full.FD.M

# combine plots
library(patchwork)
Full.FD.M/ Full.FD.T+ 
  plot_layout(heights = c(1,1), guides = "collect")->Full.FD; Full.FD

ggsave(filename = "Full FD subsets.png", plot =Full.FD, path = "results/preliminar graphs", 
       device = "png", dpi = 600)
x11()
## 5.5.1 Germ traits####
plot_x_sp_M <- as.data.frame(plot_x_sp_M)
germ.melodic.M <- melodic(plot_x_sp_M, germtraits.dist.M)
str(germ.melodic.M)

# germ melodic object data handling
as.data.frame(germ.melodic.M$abundance$mpd) %>%
  cbind (as.data.frame(germ.melodic.M$presence$mpd))%>%
  cbind (plot_x_env_M2)%>%
  mutate(Germ.Mpd.ab = germ.melodic.M$abundance$mpd)%>%
  mutate(Germ.Mpd.pa = germ.melodic.M$presence$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Germ.Mpd.ab,Germ.Mpd.pa )%>% #, Rao.pa, Mpd.pa
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, Germ.Mpd.ab,Germ.Mpd.pa , elevation, FDD, GDD, Snw)-> FD.germ.M

## 5.5.2 Plant traits Weighted  ####
weighted.plant.melodic.M <- melodic(plot_x_sp_M, Wplanttraits.dist.M)
str(weighted.plant.melodic.M)

# plant melodic object data handling
as.data.frame(weighted.plant.melodic.M$abundance$mpd) %>%
  cbind (as.data.frame(weighted.plant.melodic.M$presence$mpd))%>%
  cbind (plot_x_env_M2)%>%
  mutate(Wplant.Mpd.ab = weighted.plant.melodic.M$abundance$mpd)%>%
  mutate(Wplant.Mpd.pa = weighted.plant.melodic.M$presence$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Wplant.Mpd.ab, Wplant.Mpd.pa)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  merge(disturbance_M, by= "plot")%>%
  dplyr::select(site, plot,Wplant.Mpd.ab,Wplant.Mpd.pa,elevation, FDD, GDD, Snw, disturbance)-> FD.w.plant.M

## Join germ and plant mdp calculations
FD.germ.M %>%
  merge(FD.w.plant.M, by = c("site", "plot", "elevation", "FDD","GDD", "Snw"))-> FD.M

# first let's check MPD normality
x11()
FD.M%>%
  gather(Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = Func.div), color = "black") + 
  facet_wrap(~Func.div, ncol =2, scales = "free") +
  theme_bw(base_size = 12) # somewhat normal distributed !?!?!

## 5.5.3 MPD differences between germination and adult plant traits ####
strip_names <- c("ab"= "Abundance data", "pa"= "Presence/absence data", 
                 "Mpd" = "Mean Pairwise Dissimilarity")
FD.M%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", names = c("Traits","indices", "data_type"), too_many = "merge")%>%
  group_by(Traits, indices, data_type)%>%
  get_summary_stats(value)%>%
  mutate(Traits= as.factor(Traits))%>%
  mutate(Traits = fct_recode(Traits, "Germination"="Germ" , "Adult Plant"="Wplant"))%>%
  ggplot()+
  geom_errorbar(aes(Traits, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_point(aes(x= Traits, y= mean, fill = Traits), color = "black", shape= 21, size = 5)+
  scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  facet_grid(~data_type, labeller = as_labeller(strip_names))+
  labs(title = "Mediterranean community FD differences", y= "Mean dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")->FD.traits.M;FD.traits.M

  ggsave(FD.traits.M, file = "FD traits diff Med.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
  
FD.M %>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", 
                       names = c("Traits","indices", "data_type"), too_many = "merge")-> FD.M.diff

#PAIRED t-test?? http://www.sthda.com/english/wiki/paired-samples-t-test-in-r  

ttest.fd.diff <- function(x) {
  t.test(value ~ Traits, data = x) -> m1 #paired = TRUE, 
  broom::tidy(m1)
}

FD.M.diff%>%
  group_by(indices, data_type)%>%
  do(ttest.fd.diff(.)) %>%
  write.csv("results/FD diff MED.csv")

## 5.5.4. MPD vs microclimatic gradients#####
# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
lms.fd.env <- function(x) {
  glm(value ~ scale(elevation)+ scale(FDD) + scale(GDD) + scale(Snw) + scale (disturbance), data = x) -> m1
  broom::tidy(m1)
}
FD.M%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  group_by(Func.div)%>%
  do(lms.fd.env(.)) %>%
  write.csv("results/FD vs micro MED.csv")

ann.sig.FD.M <- data.frame (read.csv("results/sig_FD_M.csv", sep = ";"))

strip_names <- c("elevation"= "Elevation", "FDD" = "FDD", "GDD"="GDD", "Snow"="Snow")
FD.M%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(Func.div~micro_variable, scales = "free")+
  #facet_grid(factor(Func.di, levels = c())~micro_variable, scales = "free", labeller = as_labeller(strip_names))+
  geom_text(data=ann.sig.FD.M, label =ann.sig.FD.M$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Mpd in Mediterranean system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_blank(),
        legend.position = "bottom")-> FD_micro_M;FD_micro_M

ggsave(FD_micro_M, file = "FD vs micro Med.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 


## 5.5.5. explore correlation with plot richness supplementary? ####
#medFD$nbsp from scrp community level analysis
# richness x environmental gradients
as.data.frame(medFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = medFD$nbsp )%>%
  merge(FD.M, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(~micro_variable, scales = "free")+
  labs(title="Mediterranean plots richness across environmental gradients")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_env_M; richness_x_env_M

ggsave(richness_x_env_M, file = "richness x env Med.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 
# richness x FD indices
as.data.frame(medFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = medFD$nbsp )%>%
  merge(FD.M, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~Func.div, scales = "free_x", ncol=4)+
  labs(title="Mediterranean plots richness correlation with FD indices")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_FD_M; richness_x_FD_M

ggsave(richness_x_FD_M, file = "richness x FD Med.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 

# combine plots
library(patchwork)
richness_x_env_M / richness_x_FD_M + 
  plot_layout(heights = c(1,2), guides = "collect") & theme(legend.position = "bottom")-> richness_M;richness_M

ggsave(filename = "Richness_Med.png", plot =richness_M , path = "results/preliminar graphs/richness", 
       device = "png", dpi = 600)


###################################  TEMPERATE  ##################################################################
# Following chapter 5 R Trait-Based ecology by Lars####
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

# check distribution of distances germ vs plant
as.matrix(germtraits.dist.T)%>%
  as.data.frame(germtraits.dist.T)%>%
  gather(sp, germ.dist, 1:25)%>%
  cbind((as.matrix(Wplanttraits.dist.T)%>%
           as.data.frame(Wplanttraits.dist.T)%>%
           gather(sp, plant.dist, 1:25)))%>%
  dplyr::select(1,2,4)%>%
  gather(trait, dist, germ.dist:plant.dist)%>%
  filter(dist>0)%>%
  ggplot(aes(dist, fill=trait))+
  #geom_density(alpha= 0.5)+
  geom_histogram(color = "black")+
  scale_fill_manual(values = c("dodgerblue4","limegreen"))
# 5.5. Calculation of MDP (main pairwise dissimilarity) with melodic function ####
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
  rbind(MPD.T(alltraits.dist.T, plot_x_sp_T))%>%
  gather(data_type, MPD, Abundance:Presence)%>%
  mutate(trait = gsub(".dist.T", "" , trait))%>%
  group_by(trait, data_type)%>%
  get_summary_stats(MPD)%>%
  mutate(trait= as.factor(trait))%>%
  mutate(trait = fct_relevel(trait, "alltraits", "Wplanttraits", "germtraits",
                             "dark", "WP","Tconstant", 
                             "seedmass", "plantheight", "leafarea", "LDMC", "SLA"))%>%
  ggplot()+
  geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_point(aes(x= trait, y= mean, fill = trait), color = "black", shape= 21, size = 5)+
  #scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  scale_y_continuous(limits = c(0.05, 0.7))+
  coord_flip()+
  facet_grid(~data_type)+
  labs(title = "Temperate community FD differences", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        #axis.title.x = element_blank(),
        legend.position = "none")-> Full.FD.T;Full.FD.T
x11()
## 5.5.1 Germ traits####
plot_x_sp_T <- as.data.frame(plot_x_sp_T)
germ.melodic.T <- melodic(plot_x_sp_T, germ.traits.dist.T)
str(germ.melodic.T)

# germ melodic object data handling
as.data.frame(germ.melodic.T$abundance$mpd) %>%
  cbind (as.data.frame(germ.melodic.T$presence$mpd))%>%
  cbind (plot_x_env_T2)%>%
  mutate(Germ.Mpd.ab = germ.melodic.T$abundance$mpd)%>%
  mutate(Germ.Mpd.pa = germ.melodic.T$presence$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Germ.Mpd.ab, Germ.Mpd.pa)%>% #, Rao.pa, Mpd.pa
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"), by = c("plot", "elevation"))%>%
  dplyr::select(site, plot, Germ.Mpd.ab, Germ.Mpd.pa,  elevation, FDD, GDD, Snw)-> FD.germ.T

## 5.5.2 Plant traits Weighted  ####
weighted.plant.melodic.T <- melodic(plot_x_sp_T, plant.traits.weighted.dist.T)
str(weighted.plant.melodic.T)

# plant melodic object data handling
as.data.frame(weighted.plant.melodic.T$abundance$mpd) %>%
  cbind (as.data.frame(weighted.plant.melodic.T$presence$mpd))%>%
  cbind (plot_x_env_T2)%>%
  mutate(Wplant.Mpd.ab = weighted.plant.melodic.T$abundance$mpd)%>%
  mutate(Wplant.Mpd.pa = weighted.plant.melodic.T$presence$mpd)%>%
  dplyr::select(plot, elevation, FDD, GDD, Snw,Wplant.Mpd.ab,Wplant.Mpd.pa )%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"), by = c("plot", "elevation"))%>%
  merge(disturbance_T, by ="plot")%>%
  dplyr::select(site, plot,Wplant.Mpd.ab,Wplant.Mpd.pa, elevation, FDD, GDD, Snw, disturbance)-> FD.w.plant.T

## Join germ and plant mdp and rao calculations
FD.germ.T %>%
  merge(FD.w.plant.T, by = c("site", "plot", "elevation", "FDD","GDD", "Snw"))-> FD.T

# first let's check MPD and Rao normality
x11()
FD.T%>%
  gather(Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = Func.div), color = "black") + 
  facet_wrap(~Func.div, ncol =2, scales = "free") +
  theme_bw(base_size = 12) #  normally distributed 
## 5.5.3  MPD differences between germination and adult plant traits ####
strip_names <- c("ab"= "Abundance data", "pa"= "Presence/absence data", 
                 "Mpd" = "Mean Pairwise Dissimilarity")
FD.T%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", names = c("Traits","indices", "data_type"), too_many = "merge")%>%
  group_by(Traits, indices, data_type)%>%
  get_summary_stats(value)%>%
  mutate(Traits= as.factor(Traits))%>%
  mutate(Traits = fct_recode(Traits, "Germination"="Germ" , "Adult Plant"="Wplant"))%>%
  ggplot()+
  geom_errorbar(aes(Traits, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth = 0.5, size =1) +
  geom_point(aes(x= Traits, y= mean, fill = Traits), color = "black", shape= 21, size = 5)+
  scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  facet_grid(indices~data_type, labeller = as_labeller(strip_names), scales = "free_y")+
  labs(title = "Temperate community FD differences", y= "Mean dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")->FD.traits.T;FD.traits.T

ggsave(FD.traits.T, file = "FD traits diff Tem.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 

FD.T %>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", 
                       names = c("Traits","indices", "data_type"), too_many = "merge")-> FD.T.diff

ttest.fd.diff <- function(x) {
  t.test(value ~ Traits, data = x) -> m1
  broom::tidy(m1)
}

FD.T.diff%>%
  group_by(indices, data_type)%>%
  do(ttest.fd.diff(.)) %>%
  write.csv("results/FD diff TEM.csv")

## 5.5.4. MPD (weigthed by abundance) vs microclimatic gradients####

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
lms.fd.env <- function(x) {
  glm(value ~ scale(elevation) + scale(FDD) + scale(GDD) + scale(Snw), data = x) -> m1
  broom::tidy(m1)
}
FD.T%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  group_by(Func.div)%>%
  do(lms.fd.env(.)) %>%
  write.csv("results/FD vs micro TEM.csv")

ann.sig.FD.T <- data.frame (read.csv("results/sig_FD_T.csv", sep = ";"))

FD.T%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(Func.div~micro_variable, scales = "free")+
  geom_text(data=ann.sig.FD.T, label =ann.sig.FD.T$label, aes(x=x, y=y),size= 6, color = "red")+
  labs(title="Mpd in Temperate system")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_blank(),
        legend.position = "bottom")-> FD_micro_T;FD_micro_T

ggsave(FD_micro_T, file = "FD vs micro Tem.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 

## 5.5.5. explore correlation with plot richness  ####
# temFD$nbsp from scrp community level analysis
# richnes x environmental gradients
as.data.frame(temFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = temFD$nbsp )%>%
  merge(FD.T, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Rao.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(~micro_variable, scales = "free")+
  labs(title="Temperate plots richness across environmental gradients")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_env_T; richness_x_env_T

ggsave(richness_x_env_T, file = "richness x env Tem.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 
# richness x FD indices
as.data.frame(temFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = temFD$nbsp )%>%
  merge(FD.T, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Rao.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~Func.div, scales = "free_x", ncol = 4)+
  labs(title="Temperate plots richness correlation with FD indices")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_FD_T; richness_x_FD_T
                                      

ggsave(richness_x_FD_T, file = "richness x FD Tem.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 

library(patchwork)
richness_x_env_T / richness_x_FD_T + 
  plot_layout(heights = c(1,2), guides = "collect") & theme(legend.position = "bottom")-> richness_T;richness_T

ggsave(filename = "Richness_Tem.png", plot =richness_T , path = "results/preliminar graphs/richness", 
       device = "png", dpi = 600)




########## Effect size plots for functional diversity ########### #####
FD_names <- c("Germ.Mpd" = "Germ.Mpd", "Plant.Mpd" = "Plant.Mpd",
              "abundance" = "Abundace data", "presence_absence"="Presence/Absence data")


x11()
read.csv("results/lm FD estimates graphs.csv", sep=";") %>% # table with lm model results from script 7
  convert_as_factor(community, data, Func.div, term) %>%
  #mutate(trait = fct_relevel(Func.div, ))%>%
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
  facet_grid (data~Func.div,labeller = as_labeller(FD_names),  scales = "free_x") + #ncol = 8, nrow =2,
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
        axis.title.x = element_text (size=14))-> FDsig_M; FDsig_M   #FDsig_T  FDsig_M

ggsave(FDsig_M, file = "Mediterranean community FD significances.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
