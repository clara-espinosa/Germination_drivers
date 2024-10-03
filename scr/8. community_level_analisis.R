library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)

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
  mutate(trait = fct_relevel(trait, "plant_height", "seed_mass",
                             "leaf_area", "LDMC", "SLA",
                             "odds_B_dark","odds_C_WP","odds_D_constant"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CM_microclima_M
  
  for (var in unique(CM_microclima_M$micro_variable)){
    micro_plot=ggplot(CM_microclima_M[CM_microclima_M$micro_variable==var,])+
      geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
  scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~trait, ncol = 4, scales = "free")+
  labs(title="Community Means in Mediterranean system", x= var)+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")
    ggsave(micro_plot, file = paste0("CM Mediterranean x ", var,".png"),
           path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
  }

  
# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model 
summary(lm(CM_M$leaf_area~elevation+ Snw+ FDD + GDD , data=plot_x_env_M))
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cm <- function(x) {
  glm(CM ~ GDD + FDD + elevation + Snw, data = x) -> m1
  broom::tidy(m1)
}

CM_M%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_M2, by= "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms(.)) %>%
  write.csv("results/CM vs micro MED.csv")
  
# see summary in results
# we can also use the multivariate analyses for the representation, for example RDA

rda_med_all<-rda(CM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
plot(rda_med_all, type = "n", scaling = "sites")
text(rda_med_all, dis = "cn", scaling = "sites")
text(rda_med_all, dis = "sp", scaling = "sites", col = "red")

rda_med_0<-rda(CM_M~1, data= plot_x_env_M)
rda_med_all<-rda(CM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
ordistep(rda_med_0, scope=formula(rda_med_all), direction = "forward")

RsquareAdj (rda_med_all)$adj.r.squared # 0.07
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

for (var in unique(CWM_microclima_M$micro_variable)){
  micro_plot=ggplot(CWM_microclima_M[CWM_microclima_M$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Weigthed Means in Mediterranean system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CWM Mediterranean x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(CWM_M$seed_mass~elevation+ FDD + GDD + Snw, data=plot_x_env_M))
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cwm <- function(x) {
  glm(CWM ~ GDD + FDD + elevation + Snw, data = x) -> m1
  broom::tidy(m1)
}

CWM_M%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_M2, by= "plot")%>%
  gather(trait, CWM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms(.)) %>%
  write.csv("results/CWM vs micro MED.csv")

# see summary in results
# we can also use the multivariate analyses for the representation, for example RDA

rda_med_all<-rda(CWM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
plot(rda_med_all, type = "n", scaling = "sites")
text(rda_med_all, dis = "cn", scaling = "sites")
text(rda_med_all, dis = "sp", scaling = "sites", col = "red")

rda_med_0<-rda(CWM_M~1, data= plot_x_env_M)
rda_med_all<-rda(CWM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
ordistep(rda_med_0, scope=formula(rda_med_all), direction = "forward")

RsquareAdj (rda_med_all)$adj.r.squared # 0.2
# final model with elevation, FDD and Snow significant and thus included!!

### 5.3.1 CAN we trust CWM (ASK LARS to help interpretate) Spoiler NO? ####
# try only for those found significant in the LM??
# Randomization comparison ###
# consider log transform abundance data to diminish the effect of highly abundant sp
CWM_M<-functcomp(sp_x_trait_M, plot_x_sp_M, CWM.type = "all")
OBS.R2.CWM_M <- summary(lm(CWM_M$odds_C_WP ~ plot_x_env_M$Snw))$r.squared #observed R2
EXPECTED.R2.CWM_M <- vector()
for(i in 1:999) {
  random.names <- sample(rownames(sp_x_trait_M))
  sp_x_trait_M.rand1 <- sp_x_trait_M
  rownames(sp_x_trait_M.rand1) <- random.names
  sp_x_trait_M.rand.final.i <- sp_x_trait_M.rand1[sort(rownames(sp_x_trait_M.rand1)),]
  EXP.CWM_M <- functcomp(sp_x_trait_M.rand.final.i, t(sp_x_plot_M), CWM.type = "all")
  EXPECTED.R2.CWM_M[i]<-summary(lm(EXP.CWM_M$odds_C_WP ~ plot_x_env_M$Snw))$r.squared
}
head( sp_x_trait_M.rand.final.i)
hist(EXPECTED.R2.CWM_M, main = "distribution of expected R2", xlab = "")
arrows(OBS.R2.CWM_M, 250, OBS.R2.CWM_M, 40, lwd = 2)
text(OBS.R2.CWM_M, 270, "observed R2")
sum(OBS.R2.CWM_M > EXPECTED.R2.CWM_M) / 1000
# This last value represents the proportion of cases in which the observed R2 was higher than
# expected (randomly). In our case, after running the loop several times we generally found that the
# proportion was around 0.87, i.e. 87% of the cases (when considering dark and GDD)

#With this we can finally compute a p-value providing the significance of the relationship between CWM_odds_WP 
# and the GDD computed, after randomizations, as:
pval <- 1 - sum(OBS.R2.CWM_M > EXPECTED.R2.CWM_M) / 1000
pval
# In practice it tells us that the relationship between CWM_odds_dark and GDD is not significant

## 5.4 Calculation of FD (need to update with loop and functions if we are to keep this)####
# dbFD function
medFD<-dbFD(sp_x_trait_M,plot_x_sp_M, CWM.type = "all")
# we could add variable weights by the argument w, w<-c(1,1,2,2,1,3,4) for example
important.indices.M <- cbind(medFD$nbsp,medFD$FRic,medFD$FEve,medFD$FDiv,medFD$FDis,medFD$RaoQ)
colnames(important.indices.M) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
pairs(important.indices.M, pch=20)
# FRis is positively correlated with number of species as expected
# FDis and Rao very similar, actually same index with only a squating difference
# PROBLEM? Number of species positively correlated with FRic, FDis and Rao

# let's check  normality (quite normally distributed!!)
important.indices.M%>%
  as.data.frame()%>%
  gather(index, value, NumbSpecies:Rao)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = index), color = "black") + 
  facet_wrap(~index, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

data.frame(important.indices.M)%>%
  rownames_to_column(var = "plot")%>%
  gather(FDind, value, NumbSpecies:Rao)%>%
  merge(plot_x_env_M2)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(FDind = factor(FDind))%>%
  ggplot(aes(y=value, x= FDD, fill = site), color= "black")+ # GDD  FDD Snw  elevation
  geom_point(shape = 21, size =3)+
  geom_smooth (aes(y=value, x= FDD),method = "lm", se=FALSE, color = "black", inherit.aes=F)+
  facet_wrap(~FDind, nrow = 3, scales = "free")+
  labs(title="FD indices in Mediterranean system")+
  theme_classic(base_size = 12)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")

# To test the FDfor microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(medFD$FRic ~ elevation+ FDD + GDD + Snw, data=plot_x_env_M)) # nbsp   RaoQ

# 5.5 Dissimilarity matrices for MPD and RaoQ ####
# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessary with our type of data
# compute dissimilarity separately for germinations traits vs other plant traits
head(sp_x_trait_M)
sp_x_trait_M <- as.data.frame(sp_x_trait_M)
plant.traits.dist.M <- (gowdis(sp_x_trait_M["seed_mass"]) + 
                 gowdis(sp_x_trait_M["plant_height"])  +
                 gowdis(sp_x_trait_M["leaf_area"])  +
                 gowdis(sp_x_trait_M["LDMC"])  +
                 gowdis(sp_x_trait_M["SLA"])) / 5

germ.traits.dist.M <- (gowdis(sp_x_trait_M["odds_B_dark"]) +
                          gowdis(sp_x_trait_M["odds_C_WP"]) +
                          gowdis(sp_x_trait_M["odds_D_constant"])) / 3

FD.planttraits.M <- dbFD(plant.traits.dist.M, plot_x_sp_M, message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.planttraits.M$RaoQ ~ elevation + FDD + GDD + Snw, data=plot_x_env_M))
# GDD and Snow significant
FD.germtraits.M <- dbFD(germ.traits.dist.M, plot_x_sp_M, message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.germtraits.M$RaoQ ~ elevation + FDD + GDD + Snw, data=plot_x_env_M))
# only elevation significantly affect germination RaoQ

### 5.5.1 cluster dendograma (does not make a lot of sense, remove??) ####
tree.traits.M <- hclust(germ.traits.dist.M, "average")
plot(tree.traits.M)

ultra.tree.traits.M <- as.phylo(tree.traits.M)
plot(ultra.tree.traits.M)

faith.FD.M <- pd(plot_x_sp_M,ultra.tree.traits.M )
FD.alltraits.M <- dbFD(all.dist.M, log(plot_x_sp_M+1), message = F, calc.CWM =F, stand.FRic = T)

par(mfrow=c(1,2))
par(mar = c(4, 4, 2, 1))
plot(faith.FD.M$SR, faith.FD.M$PD, pch=20, ylab="Faith FD", xlab="Number of species")
plot(FD.alltraits.M$FRic, faith.FD.M$PD, pch=20, ylab="Faith FD", xlab="Functional Richness")
# these indicies provide redundant information with the number of species
# to remove the effect of species richness on FD indicies use function ses.pd in picante package
# by default the function will randomize species name and use 999 randomizations
ses.faith.FD.M <- ses.pd(plot_x_sp_M, ultra.tree.traits.M, runs = 199)

plot(ses.faith.FD.M$ntaxa, ses.faith.FD.M$pd.obs.z, pch = 20, # ntaxa = speceis richness
     ylab = "SES Faith FD", xlab = "Number of species")
plot(FD.alltraits.M$RaoQ, ses.faith.FD.M$pd.obs.z, pch = 20,  # pd.obs.z = standarized values
     ylab = "SES Faith FD", xlab = "Rao")
par(mfrow=c(1,1))
# FAITH index no longer correlated to number of species, but it generally reflects quite well the Rao index 

## 5.6 START FROM HERE!! Calculation of MDP (main pairwise dissimilarity) and RAo with picante package ####
## 5.6.1 Germination traits (there is an NA, why? all traits are complete)####
germ.mpd.M <- mpd(plot_x_sp_M, as.matrix(germ.traits.dist.M))
## Species abundances can be considered with the argument abundance.weighted = TRUE
# we could logt ransform abundance data to damp a bit the effect of overabundant species
germ.mpd.M.ab <- mpd(plot_x_sp_M, as.matrix(germ.traits.dist.M), abundance.weighted = T)

plot_x_sp_M <- data.frame(plot_x_sp_M)
Rao.dvc.M <- divc(plot_x_sp_M, germ.traits.dist.M)
# error in divc, non coveniant df

# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
#Error in if (class(dis) != "matrix") { : the condition has length > 1????
plot_x_sp_M <- as.matrix(plot_x_sp_M)
germ.melodic.M <- melodic(plot_x_sp_M, as.matrix(germ.traits.dist.M))
str(germ.melodic.M)
# MPD = Rao/Simpson

## 5.6.2 Other plant traits ####
plant.mpd.M <- mpd(plot_x_sp_M, as.matrix(plant.traits.dist.M))
## Species abundances can be considered with the argument abundance.weighted = TRUE
# we could log transform abundance data to damp a bit the effect of overabundant species
plant.mpd.M.ab <- mpd(plot_x_sp_M, as.matrix(plant.traits.dist.M), abundance.weighted = T)

plot_x_sp_M <- data.frame(plot_x_sp_M)
Rao.dvc.M <- divc(plot_x_sp_M, plant.traits.dist.M)
# error in divc, non coveniant df

# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
#Error in if (class(dis) != "matrix") { : the condition has length > 1????
plot_x_sp_M <- as.matrix(plot_x_sp_M)
plant.melodic.M <- melodic(plot_x_sp_M, as.matrix(plant.traits.dist.M))
str(plant.melodic.M)
# MPD = Rao/Simpso

## below not used ####
mntd.M <- mntd(plot_x_sp_M, as.matrix(all.dist.M))
number.species.M <- medFD$nbsp

par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
plot(number.species.M, mpd.M, pch = 20, ylab = "MPD",
     xlab = "Number of species")
plot(number.species.M, mntd.M, pch = 20, ylab = "MNTD",
     xlab = "Number of species")
par(mfrow = c(1, 1))
# We observe a general lack of correlation between MPD and the number of species
# as expected, we observe a negative correlation between MNTD and the number of species

## Species abundances can be considered with the argument abundance.weighted = TRUE
# we could logtransform abundance data to damp a bit the effect of overabundant species

mpd.M.ab <- mpd(plot_x_sp_M, as.matrix(all.dist.M), abundance.weighted = T)
mntd.M.ab <- mntd(plot_x_sp_M, as.matrix(all.dist.M), abundance.weighted = T)

plot(FD.alltraits.M$RaoQ, mpd.M.ab, pch=20, ylab= "MPd abundance")

# not working
plot_x_sp_M <- data.frame(plot_x_sp_M)
Rao.dvc.M <- divc(plot_x_sp_M, all.dist.M)

# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
#Error in if (class(dis) != "matrix") { : the condition has length > 1????
plot_x_sp_M <- as.matrix(plot_x_sp_M)
melodic.M <- melodic(plot_x_sp_M, as.matrix(all.dist.M))
str(melodic.M)
# MPD = Rao/Simpson

plot(melodic.M$abundance$mpd, melodic.M$abundance$rao, 
     xlab="MPD", ylab="Rao", pch=20, main= "Mediterranean community (species poor)")
abline(0, 1)

# Use of Rao function in scr folder ####
Rao.M <- Rao(sp_x_plot_M, all.dist.M, dphyl = NULL, weight = F, Jost=T, structure = NULL)

Rao.M$TD$Mean_Alpha  # alpha taxonomic diversity 
Rao.M$TD$Beta_add  # beta taxonomic diversity (additive) 
Rao.M$TD$Beta_prop  # beta taxonomic diversity (proportional) 
Rao.M$TD$Gamma  # alpha taxonomic diversity 

round(Rao.M$TD$Pairwise_samples$Beta_prop) # without standarization
round(Rao.M$TD$Pairwise_samples$Beta_prop/(1-1/73)) # 73 is the total of plots?

par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
TD.M <- c(Rao.M$TD$Mean_Alpha, Rao.M$TD$Beta_add)
FD.M <- c(Rao.M$FD$Mean_Alpha, Rao.M$FD$Beta_add)
barplot(cbind(TD, FD.M), beside = F, col = c(8, 0),
        ylim = c(0, 7), ylab = c("Rao Equivalent Numbers"), main= "Mediterranean community")
points(2, 4, pch = 22, bg = "grey")
points(2, 3.5, pch = 22, bg = "white")
text(2, 4, "alpha", pos = 4)
text(2, 3.5, "beta", pos = 4)

TD.Rao.M.prop <- c(100 * Rao.M$TD$Mean_Alpha / Rao.M$TD$Gamma,
                   Rao.M$TD$Beta_prop)
FD.Rao.M.prop <- c(100*Rao.M$FD$Mean_Alpha/Rao.M$FD$Gamma,
                   Rao.M$FD$Beta_prop)
barplot(cbind(TD.Rao.M.prop, FD.Rao.M.prop), beside = F, col = c(8, 0),
        ylim = c(0, 100), ylab = c("Proportion of Alpha and Beta in Equivalent numbers"))
text(0.7, 60, "beta")
text(0.7, 15, "alpha")
par(mfrow = c(1, 1))
#despite the great turnover in species composition (beta TD, grey colour, left stack) 
# we have a real low functional turnover. In other words most of the trait dissimilarity between 
# species is found within a plot and not across plots, even if there are such a marked 
# environmental changes and even if species composition changes are very strong (high beta TD)

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
# 5.2 Calculation of CM ####
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

# loop for visualization
CM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "plant_height", "seed_mass", 
                             "leaf_area", "LDMC", "SLA",
                             "odds_B_dark","odds_C_WP","odds_D_constant"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CM_microclima_T

for (var in unique(CM_microclima_T$micro_variable)){
  micro_plot=ggplot(CM_microclima_T[CM_microclima_T$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Means in Temperate system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CM Temperate x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}


# GDD, FDD, Snow, Elevation
# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(CM_T$seed_production~elevation+ Snw+ FDD + GDD , data=plot_x_env_T))
# see summary in results
# (function written to repeat each lm x microclimatic variable)
# formula adapted from PICOS github repository scr Fig5-GLMs
lms.cm <- function(x) {
  glm(CM ~ GDD + FDD + elevation + Snw, data = x) -> m1
  broom::tidy(m1)
}

CM_T%>%
  rownames_to_column(var = "plot")%>%
  merge(plot_x_env_T2, by= "plot")%>%
  gather(trait, CM, seed_mass:odds_D_constant)%>%
  group_by (trait)%>%
  do(lms.cm(.)) %>%
  write.csv("results/CM vs micro TEM.csv")
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

# loop for visualization
CWM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:odds_D_constant)%>%
  merge(plot_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(trait = factor(trait))%>%
  mutate(trait = fct_relevel(trait, "plant_height", "seed_mass", 
                             "leaf_area", "LDMC", "SLA",
                             "odds_B_dark","odds_C_WP","odds_D_constant"))%>%
  dplyr::select(site, plot, trait, value, elevation, Snw, FDD, GDD)%>%
  rename(Snow=Snw)%>%
  gather(micro_variable, value_micro,  elevation:GDD)-> CWM_microclima_T

for (var in unique(CWM_microclima_T$micro_variable)){
  micro_plot=ggplot(CWM_microclima_T[CWM_microclima_T$micro_variable==var,])+
    geom_point(aes(y=value, x= value_micro, fill = site), color= "black",shape = 21, size =4)+
    scale_fill_manual( values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
    geom_smooth (aes(y=value, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
    facet_wrap(~trait, ncol = 4, scales = "free")+
    labs(title="Community Weighted Means in Temperate system", x= var)+
    theme_classic(base_size = 16)+
    theme(strip.text = element_text(size =12),
          legend.position = "bottom")
  ggsave(micro_plot, file = paste0("CWM Temperate x ", var,".png"),
         path = "results/preliminar graphs", scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(CWM_T$seed_production~elevation+ FDD + GDD + Snw, data=plot_x_env_T))

lms.cwm <- function(x) {
  glm(CWM ~ GDD + FDD + elevation + Snw, data = x) -> m1
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
# we can also use the multivariate analyses for the representation, for example RDA

rda_tem_all<-rda(CWM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
 plot(rda_tem_all, type = "n", scaling = "sites")
 text(rda_tem_all, dis = "cn", scaling = "sites")
 text(rda_tem_all, dis = "sp", scaling = "sites", col = "red")

rda_tem_0<-rda(CWM_T~1, data= plot_x_env_T)
rda_tem_all<-rda(CWM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
ordistep(rda_tem_0, scope=formula(rda_tem_all), direction = "forward")

RsquareAdj (rda_tem_all)$adj.r.squared # 0.03
 # rda_tem_0 is the best model!

## 5.3.1 CAN we trust CWM (spoiler NO) ?? ####
# Randomization comparison ###
CWM_T<-functcomp(sp_x_trait_T, t(sp_x_plot_T), CWM.type = "all")
OBS.R2.CWM_T <- summary(lm(CWM_T$plant_height ~ plot_x_env_T$GDD))$r.squared #observed R2
EXPECTED.R2.CWM_T <- vector()
for(i in 1:999) {
  random.names <- sample(rownames(sp_x_trait_T))
  sp_x_trait_T.rand1 <- sp_x_trait_T
  rownames(sp_x_trait_T.rand1) <- random.names
  sp_x_trait_T.rand.final.i <- sp_x_trait_T.rand1[sort(rownames(sp_x_trait_T.rand1)),]
  EXP.CWM_T <- functcomp(sp_x_trait_T.rand.final.i, t(sp_x_plot_T), CWM.type = "all")
  EXPECTED.R2.CWM_T[i]<-summary(lm(EXP.CWM_T$plant_height ~ plot_x_env_T$GDD))$r.squared
}
head( sp_x_trait_T.rand.final.i)
hist(EXPECTED.R2.CWM_T, main = "distribution of expected R2", xlab = "")
arrows(OBS.R2.CWM_T, 250, OBS.R2.CWM_T, 40, lwd = 2)
text(OBS.R2.CWM_T, 270, "observed R2")
sum(OBS.R2.CWM_T > EXPECTED.R2.CWM_T) / 1000
# This last value represents the proportion of cases in which the observed R2 was higher than
# expected. In our case, after running the loop several times we generally found that the
# proportion was around 0.79, i.e. 79% of the cases

#With this we can finally compute a p-value providing the significance of the relationship between CWM_odds_WP 
# and the GDD computed, after randomizations, as:
pval <- 1 - sum(OBS.R2.CWM_T > EXPECTED.R2.CWM_T) / 1000
pval
# In practice it tells us that the relationship between CWM_odds_WP and elevation is not significant

# 5.4 Calculation of FD ####
# dbFD function
temFD<-dbFD(sp_x_trait_T,plot_x_sp_T, CWM.type = "all")
# we could add variable weights by the argument w, w<-c(1,1,2,2,1,3,4) for example
important.indices.T <- cbind(temFD$nbsp,temFD$FRic,temFD$FEve,temFD$FDiv,temFD$FDis,temFD$RaoQ)
colnames(important.indices.T) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
pairs(important.indices.T, pch=20)
# FRis is positively correlated with number of species as expected
# FDis and Rao very similar, actually same index with only a squating difference

# let's check  normality (quite normally distributed!!)
important.indices.T%>%
  as.data.frame()%>%
  gather(index, value, NumbSpecies:Rao)%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = index), color = "black") + 
  facet_wrap(~index, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

data.frame(important.indices.T)%>%
  rownames_to_column(var = "plot")%>%
  gather(FDind, value, NumbSpecies:Rao)%>%
  merge(plot_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(FDind = factor(FDind))%>%
  ggplot(aes(y=value, x= Snw, fill = site), color= "black")+ # GDD  FDD Snw  elevation
  geom_point(shape = 21, size =3)+
  geom_smooth (aes(y=value, x= Snw),method = "lm", se=FALSE, color = "black", inherit.aes=F)+
  facet_wrap(~FDind, nrow = 3, scales = "free")+
  labs(title="FD indices in Temperate system")+
  theme_classic(base_size = 12)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(temFD$FDis~elevation+ FDD + GDD + Snw, data=plot_x_env_T)) # nbsp   RaoQ

# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessaty with our type of data
head(sp_x_trait_T)
sp_x_trait_T <- as.data.frame(sp_x_trait_T)
all.dist.T <- (gowdis(sp_x_trait_T["seed_mass"]) + gowdis(sp_x_trait_T["seed_production"]) + 
                 gowdis(sp_x_trait_T["plant_height"]) + gowdis(sp_x_trait_T["floral_height"]) +
                 gowdis(sp_x_trait_T["A_control"]) +gowdis(sp_x_trait_T["B_dark"]) +
                 gowdis(sp_x_trait_T["C_WP"]) +gowdis(sp_x_trait_T["D_constant"]) +
                 gowdis(sp_x_trait_T["odds_B_dark"]) + gowdis(sp_x_trait_T["E_cold"])+
                 gowdis(sp_x_trait_T["odds_C_WP"]) +gowdis(sp_x_trait_T["odds_D_constant"]) +
                 gowdis(sp_x_trait_T["odds_E_cold"])) / 9
FD.alltraits.T <- dbFD(all.dist.T, log(plot_x_sp_T+1), message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.alltraits.T$RaoQ ~elevation+ FDD + GDD + Snw, data=plot_x_env_T))

# 5.5 Dissimilarity matrices for MPD and RaoQ ####
# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessary with our type of data
# compute dissimilarity separately for germinations traits vs other plant traits
head(sp_x_trait_T)
sp_x_trait_T <- as.data.frame(sp_x_trait_T)
plant.traits.dist.T <- (gowdis(sp_x_trait_T["seed_mass"]) + 
                          gowdis(sp_x_trait_T["plant_height"])  +
                          gowdis(sp_x_trait_T["leaf_area"])  +
                          gowdis(sp_x_trait_T["LDMC"])  +
                          gowdis(sp_x_trait_T["SLA"])) / 5

germ.traits.dist.T <- (gowdis(sp_x_trait_T["odds_B_dark"]) +
                         gowdis(sp_x_trait_T["odds_C_WP"]) +
                         gowdis(sp_x_trait_T["odds_D_constant"])) / 3

FD.planttraits.T <- dbFD(plant.traits.dist.T, plot_x_sp_T, message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.planttraits.T$RaoQ ~ elevation + FDD + GDD + Snw, data=plot_x_env_T))
# GDD significant
FD.germtraits.T <- dbFD(germ.traits.dist.T, plot_x_sp_T, message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.germtraits.T$RaoQ ~ elevation + FDD + GDD + Snw, data=plot_x_env_T))
# nothing significant
## 5.5.1 cluster dendograma (does not make a lot of sense) ####
tree.traits.T <- hclust(germ.traits.dist.T, "average")
plot(tree.traits.T)

ultra.tree.traits.T <- as.phylo(tree.traits.T)
plot(ultra.tree.traits.T)

faith.FD.T <- pd(plot_x_sp_T,ultra.tree.traits.T )
FD.alltraits.T <- dbFD(all.dist.T, log(plot_x_sp_T+1), message = F, calc.CWM =F, stand.FRic = T)

par(mfrow=c(1,2))
par(mar = c(4, 4, 2, 1))
plot(faith.FD.T$SR, faith.FD.T$PD, pch=20, ylab="Faith FD", xlab="Number of species")
plot(FD.alltraits.T$FRic, faith.FD.T$PD, pch=20, ylab="Faith FD", xlab="Functional Richness")
# these indicies provide redundant information with the number of species
# to remove the effect of species richness on FD indicies use function ses.pd in picante package
# by default the function will randomize species name and use 999 randomizations
ses.faith.FD.T <- ses.pd(plot_x_sp_T, ultra.tree.traits.T, runs = 199)
par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
plot(ses.faith.FD.T$ntaxa, ses.faith.FD.T$pd.obs.z, pch = 20, # ntaxa = speceis richness
     ylab = "SES Faith FD", xlab = "Number of species")
plot(FD.alltraits.T$RaoQ, ses.faith.FD.T$pd.obs.z, pch = 20,  # pd.obs.z = standarized values
     ylab = "SES Faith FD", xlab = "Rao")
# FAITH index no longer correlated to number of species, but it generally reflects quite well the Rao index 

## 5.6 START FROM HERE!! Calculation of MDP (main pairwise dissimilarity) and RAo with picante package ####
mpd.T <- mpd(plot_x_sp_T, as.matrix(all.dist.T))
mntd.T <- mntd(plot_x_sp_T, as.matrix(all.dist.T))
number.species.T <- temFD$nbsp

par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
plot(number.species.T, mpd.T, pch = 20, ylab = "MPD",
     xlab = "Number of species")
plot(number.species.T, mntd.T, pch = 20, ylab = "MNTD",
     xlab = "Number of species")
par(mfrow = c(1, 1))
# We observe a general lack of correlation between MPD and the number of species
# as expected, we observe a negative correlation between MNTD and the number of species
## Species abundances can be considered with the argument abundance.weighted = TRUE
# we could logtransforme abundance data to damp a bit the effect of overabundant species

mpd.T.ab <- mpd(plot_x_sp_T, as.matrix(all.dist.T), abundance.weighted = T)
mntd.T.ab <- mntd(plot_x_sp_T, as.matrix(all.dist.T), abundance.weighted = T)

plot(FD.alltraits.T$RaoQ, mpd.T.ab, pch=20, ylab= "MPd abundance")

# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
# Error in if (class(samp) != "matrix") { : the condition has length > 1
melodic.T <- melodic(plot_x_sp_T, as.matrix(all.dist.T))
str(melodic.T)

# MPD = Rao/Simpson
par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
plot(melodic.T$abundance$mpd, melodic.T$abundance$rao, 
     xlab="MPD", ylab="Rao", pch=20, main= "Temperate community (species rich)")
abline(0, 1)

# Use of Rao function in scr folder
Rao.T <- Rao(sp_x_plot_T, all.dist.T, dphyl = NULL, weight = F, Jost=T, structure = NULL)

Rao.T$TD$Mean_Alpha  # alpha taxonomic diversity 
Rao.T$TD$Beta_add  # beta taxonomic diversity (additive) 
Rao.T$TD$Beta_prop  # beta taxonomic diversity (proportional) 
Rao.T$TD$Gamma  # alpha taxonomic diversity 

round(Rao.T$TD$Pairwise_samples$Beta_prop) # without standarization
round(Rao.T$TD$Pairwise_samples$Beta_prop/(1-1/73)) # 73 is the total of plots?

TD.T <- c(Rao.T$TD$Mean_Alpha, Rao.T$TD$Beta_add)
FD.T <- c(Rao.T$FD$Mean_Alpha, Rao.T$FD$Beta_add)
barplot(cbind(TD.T, FD.T), beside = F, col = c(8, 0),
        ylim = c(0, 17), ylab = c("Rao Equivalent Numbers"), main= "Temperate community")
points(2, 4, pch = 22, bg = "grey")
points(2, 3.5, pch = 22, bg = "white")
text(2, 4, "alpha", pos = 4)
text(2, 3.5, "beta", pos = 4)

TD.Rao.T.prop <- c(100 * Rao.T$TD$Mean_Alpha / Rao.T$TD$Gamma,
                   Rao.T$TD$Beta_prop)
FD.Rao.T.prop <- c(100*Rao.T$FD$Mean_Alpha/Rao.T$FD$Gamma,
                   Rao.T$FD$Beta_prop)
barplot(cbind(TD.Rao.T.prop, FD.Rao.T.prop), beside = F, col = c(8, 0),
        ylim = c(0, 100), ylab = c("Proportion of Alpha and Beta in Equivalent numbers"))
text(0.7, 60, "beta")
text(0.7, 15, "alpha")
#despite the great turnover in species composition (beta TD, grey colour, left stack) 
# we have a real low functional turnover. In other words most of the trait dissimilarity between 
# species is found within a plot and not across plots, even if there are such a marked 
# environmental changes and even if species composition changes are very strong (high beta TD)
