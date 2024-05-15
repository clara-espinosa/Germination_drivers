library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
Sys.setlocale("LC_TIME", "English")

# preliminary steps #
###################################  MEDITERRANEAN ###############################################################
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

setdiff(spatial_sp_Med$species, species$species)
# 2A- Cover/plot of species with germination traits ####
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

# 2B- Number of species with germination traits per plot ####
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

# 3- temperature graphs of plots with 80% coverage with trait data ####
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

# CWM test  Following chapter 5 R TRait-Based ecology by Lars####
# 5.1 Data matrices ####
# in sp x plot 
# in all files names of species and plots are not part of the matrix, used as row and columns names
# mediterranean #
sp_x_plot_M # species in rows and plots in columns, abundance already as relative cover
plot_x_sp_M  # plot in row and species in columns, abundance in relative cover
sp_x_trait_M # species in rows and traits in columns
# Consider remove species which had 0 germination across al treatments in germination experiment
# Solidago virgaurea
plot_x_env_M # plot in rows and env data in columns

## 5.2 Check normality of the quantitative traits ####
### all traits are right skewed, try a log transformation look better but still not normal
sp_x_trait_M%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

## 5.3 Calculation of CWM ####
# consider log transform abundances to reduce affects of very abundant species
CWM_M<-functcomp(sp_x_trait_M, plot_x_sp_M, CWM.type = "all") # CWM.type to consider binary or categorical traits
CWM_M<-functcomp(sp_x_trait_M, log(plot_x_sp_M+1), CWM.type = "all")
# Let’s see if CWM changes along the climatic gradient
# explore visually a bit the results before running any analysis,
plot_x_env_M%>%
  rownames_to_column(var = "plot")->plot_x_env_M2 
# germination proportion traits probably shouldn't be considered 
CWM_M%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  merge(plot_x_env_M2)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(trait = factor(trait))%>%
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
# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(CWM_M$E_cold_stratification~elevation+ FDD + GDD + Snw, data=plot_x_env_M))
# see summary in results
# we can also use the multivariate analyses for the representation, for example RDA

rda_med_all<-rda(CWM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
plot(rda_med_all, type = "n", scaling = "sites")
text(rda_med_all, dis = "cn", scaling = "sites")
text(rda_med_all, dis = "sp", scaling = "sites", col = "red")

rda_med_0<-rda(CWM_M~1, data= plot_x_env_M)
rda_med_all<-rda(CWM_M~elevation+ FDD + GDD + Snw , data= plot_x_env_M)
ordistep(rda_med_0, scope=formula(rda_med_all), direction = "forward")

RsquareAdj (rda_med_all)$adj.r.squared # 0.06
# final model only GDD is significant and thus included!!

### 5.3.1 CAN we trust CWM ?? ####
# Randomization comparison ###
# consider log transform abundance data to diminish the effect of highly abundant sp
CWM_M<-functcomp(sp_x_trait_M, plot_x_sp_M, CWM.type = "all")
OBS.R2.CWM_M <- summary(lm(CWM_M$odds_C_WP ~ plot_x_env_M$GDD))$r.squared #observed R2
EXPECTED.R2.CWM_M <- vector()
for(i in 1:999) {
  random.names <- sample(rownames(sp_x_trait_M))
  sp_x_trait_M.rand1 <- sp_x_trait_M
  rownames(sp_x_trait_M.rand1) <- random.names
  sp_x_trait_M.rand.final.i <- sp_x_trait_M.rand1[sort(rownames(sp_x_trait_M.rand1)),]
  EXP.CWM_M <- functcomp(sp_x_trait_M.rand.final.i, t(sp_x_plot_M), CWM.type = "all")
  EXPECTED.R2.CWM_M[i]<-summary(lm(EXP.CWM_M$odds_C_WP ~ plot_x_env_M$GDD))$r.squared
}
head( sp_x_trait_M.rand.final.i)
hist(EXPECTED.R2.CWM_M, main = "distribution of expected R2", xlab = "")
arrows(OBS.R2.CWM_M, 250, OBS.R2.CWM_M, 40, lwd = 2)
text(OBS.R2.CWM_M, 270, "observed R2")
sum(OBS.R2.CWM_M > EXPECTED.R2.CWM_M) / 1000
# This last value represents the proportion of cases in which the observed R2 was higher than
# expected. In our case, after running the loop several times we generally found that the
# proportion was around 0.64, i.e. 64% of the cases

#With this we can finally compute a p-value providing the significance of the relationship between CWM_odds_WP 
# and the GDD computed, after randomizations, as:
pval <- 1 - sum(OBS.R2.CWM_M > EXPECTED.R2.CWM_M) / 1000
pval
# In practice it tells us that the relationship between CWM_odds_WP and GDD is not significant

## 5.4 Calculation of FD ####
# dbFD function
medFD<-dbFD(sp_x_trait_M,log(plot_x_sp_M+1), CWM.type = "all")
# we could add variable weights by the argument w, w<-c(1,1,2,2,1,3,4) for example
important.indices.M <- cbind(medFD$nbsp,medFD$FRic,medFD$FEve,medFD$FDiv,medFD$FDis,medFD$RaoQ)
colnames(important.indices.M) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
pairs(important.indices.M, pch=20)
# FRis is positively correlated with number of species as expected
# FDis and Rao very similar, actually same index with only a squating difference
# PROBLEM? Number of species positively correlated with FRic, FDis and Rao
data.frame(important.indices.M)%>%
  rownames_to_column(var = "plot")%>%
  gather(FDind, value, NumbSpecies:Rao)%>%
  merge(plot_x_env_M2)%>%
  merge(read.csv("data/spatial-survey-header-Med.csv"), by = c("plot", "elevation"))%>%
  mutate(FDind = factor(FDind))%>%
  ggplot(aes(y=value, x= Snw, fill = site), color= "black")+ # GDD  FDD Snw  elevation
  geom_point(shape = 21, size =3)+
  geom_smooth (aes(y=value, x= Snw),method = "lm", se=FALSE, color = "black", inherit.aes=F)+
  facet_wrap(~FDind, nrow = 3, scales = "free")+
  labs(title="FD indices in Mediterranean system")+
  theme_classic(base_size = 12)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(medFD$RaoQ ~ elevation+ FDD + GDD + Snw, data=plot_x_env_M)) # nbsp   RaoQ

# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessaty with our type of data
head(sp_x_trait_M)
all.dist.M <- (gowdis(sp_x_trait_M["seed_mass"]) + gowdis(sp_x_trait_M["seed_production"]) + 
                 gowdis(sp_x_trait_M["plant_height"]) + gowdis(sp_x_trait_M["floral_height"]) +
                 gowdis(sp_x_trait_M["odds_A_control"]) +gowdis(sp_x_trait_M["odds_B_dark"]) +
                 gowdis(sp_x_trait_M["odds_C_WP"]) +gowdis(sp_x_trait_M["odds_D_constant"]) +
                 gowdis(sp_x_trait_M["odds_E_cold"])) / 9
FD.alltraits.M <- dbFD(all.dist.M, log(plot_x_sp_M+1), message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.alltraits.M$RaoQ ~elevation+ FDD + GDD + Snw, data=plot_x_env_M))

### 5.4.1 cluster dendograma ####
tree.traits.M <- hclust(all.dist.M, "average")
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

## 5.5 Calculation of MDP (main pairwise dissimilarity) and RAo with picante package ####
mpd.M <- mpd(plot_x_sp_M, as.matrix(all.dist.M))
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
# we could logtransforme abundance data to damp a bit the effect of overabundant species

mpd.M.ab <- mpd(log(plot_x_sp_M+1), as.matrix(all.dist.M), abundance.weighted = T)
mntd.M.ab <- mntd(log(plot_x_sp_M+1), as.matrix(all.dist.M), abundance.weighted = T)

plot(FD.alltraits.M$RaoQ, mpd.M.ab, pch=20, ylab= "MPd abundance")

# not working
Rao.dvc.M <- divc(plot_x_sp_M, all.dist.M)

# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
melodic.M <- melodic(log(plot_x_sp_M+1), as.matrix(all.dist.M))
str(melodic.M)
# MPD = Rao/Simpson

plot(melodic.M$abundance$mpd, melodic.M$abundance$rao, 
     xlab="MPD", ylab="Rao", pch=20, main= "Mediterranean community (species poor)")
abline(0, 1)

# Use of Rao function in scr folder
Rao.M <- Rao(log(sp_x_plot_M+1), all.dist.M, dphyl = NULL, weight = F, Jost=T, structure = NULL)

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
# 1-check species names to match####
# dataframe with species data from germination experiments
read.csv("data/species.csv", sep=",") %>%
  select(species, code,  family, community, habitat, germ_drivers)  %>%
  convert_as_factor(species, code, family, community, habitat, germ_drivers)-> species 

# dataframe with species data from spatial survay
read.csv("data/spatial-survey-species-Tem.csv", sep=",")%>%
  group_by(plot)%>%
  mutate(total_cover= sum(cover),
         rel_cover = (cover/total_cover)*100,
         N_sp = length(unique(species)))-> spatial_sp_Tem
setdiff(spatial_sp_Tem$species, species$species)

# 2A- Cover/plot of species with germination traits ####
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
species %>%
  #filter(germ_drivers == "Yes")%>%
  merge(spatial_sp_Tem, by = "species")%>%
  #filter (!species == "Euphrasia salisburgensis")%>% # species with 0 germ across all treatments
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
tem_plot1%>%
  merge(read.csv("data/spatial-survey-temperatures-Tem.csv")) %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC", format = "%d/%m/%Y %H:%M" )) %>% 
  merge(read.csv("data/spatial-survey-header-Tem.csv"))-> spatial1 

setdiff(tem_plot1$plot, spatial1$plot)
unique(tem_plot1$plot)
unique(spatial1$plot)

spatial1 %>% 
  group_by(Time = lubridate::floor_date(Time, "day"), site, plot) %>%
  summarise(Temperature = mean(Temperature)) -> spatial2

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

sp_x_plot_T# species in rows and plots in columns, abundance already as relative cover
plot_x_sp_T# plot in row and species in columns, abundance in relative cover
sp_x_trait_T# species in rows and traits in columns
# Consider remove species which had 0 germination across al treatments in germination experiment
plot_x_env_T# plot in rows and env data in columns

## 5.2 Check normality of the quantitative traits ####
### all traits are right skewed, ty a log transformation look better but still not normal
sp_x_trait_T%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  #mutate(value = log(value))%>%
  ggplot()+
  geom_histogram(aes(x= value, fill = trait), color = "black") + 
  facet_wrap(~trait, ncol = 3, scales = "free") +
  theme_bw(base_size = 12)

## 5.3 Calculation of CWM ####
# consider log transform abundance data 
CWM_T<-functcomp(sp_x_trait_T, t(sp_x_plot_T), CWM.type = "all") # CWM.type to consider binary or categorical traits
CWM_T<-functcomp(sp_x_trait_T, plot_x_sp_T, CWM.type = "all") 
# Let’s see if CWM changes along the climatic gradient
# explore visually a bit the results before running any analysis,
plot_x_env_T%>%
  rownames_to_column(var = "plot")->plot_x_env_T2 

CWM_T%>%
  rownames_to_column(var = "plot")%>%
  gather(trait, value, seed_mass:E_cold_stratification)%>%
  merge(plot_x_env_T2)%>%
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
summary(lm(CWM_T$E_cold_stratification~elevation+ FDD + GDD + Snw, data=plot_x_env_T))
# see summary in results
# we can also use the multivariate analyses for the representation, for example RDA

rda_tem_all<-rda(CWM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
 plot(rda_tem_all, type = "n", scaling = "sites")
 text(rda_tem_all, dis = "cn", scaling = "sites")
 text(rda_tem_all, dis = "sp", scaling = "sites", col = "red")

rda_tem_0<-rda(CWM_T~1, data= plot_x_env_T)
rda_tem_all<-rda(CWM_T~elevation+ FDD + GDD + Snw , data= plot_x_env_T)
ordistep(rda_tem_0, scope=formula(rda_tem_all), direction = "forward")
 # final model only elevation is significant and thus included!!
## 5.3.1 CAN we trust CWM ?? ####
# Randomization comparison ###
CWM_T<-functcomp(sp_x_trait_T, t(sp_x_plot_T), CWM.type = "all")
OBS.R2.CWM_T <- summary(lm(CWM_T$odds_C_WP ~ plot_x_env_T$elevation))$r.squared #observed R2
EXPECTED.R2.CWM_T <- vector()
for(i in 1:999) {
  random.names <- sample(rownames(sp_x_trait_T))
  sp_x_trait_T.rand1 <- sp_x_trait_T
  rownames(sp_x_trait_T.rand1) <- random.names
  sp_x_trait_T.rand.final.i <- sp_x_trait_T.rand1[sort(rownames(sp_x_trait_T.rand1)),]
  EXP.CWM_T <- functcomp(sp_x_trait_T.rand.final.i, t(sp_x_plot_T), CWM.type = "all")
  EXPECTED.R2.CWM_T[i]<-summary(lm(EXP.CWM_T$odds_C_WP ~ plot_x_env_T$elevation))$r.squared
}
head( sp_x_trait_T.rand.final.i)
hist(EXPECTED.R2.CWM_T, main = "distribution of expected R2", xlab = "")
arrows(OBS.R2.CWM_T, 250, OBS.R2.CWM_T, 40, lwd = 2)
text(OBS.R2.CWM_T, 270, "observed R2")
sum(OBS.R2.CWM_T > EXPECTED.R2.CWM_T) / 1000
# This last value represents the proportion of cases in which the observed R2 was higher than
# expected. In our case, after running the loop several times we generally found that the
# proportion was around 0.9, i.e. 90% of the cases

#With this we can finally compute a p-value providing the significance of the relationship between CWM_odds_WP 
# and the GDD computed, after randomizations, as:
pval <- 1 - sum(OBS.R2.CWM_T > EXPECTED.R2.CWM_T) / 1000
pval
# In practice it tells us that the relationship between CWM_odds_WP and elevation is not significant

## 5.4 Calculation of FD ####
# dbFD function
temFD<-dbFD(sp_x_trait_T,log(plot_x_sp_T+1), CWM.type = "all")
# we could add variable weights by the argument w, w<-c(1,1,2,2,1,3,4) for example
important.indices.T <- cbind(temFD$nbsp,temFD$FRic,temFD$FEve,temFD$FDiv,temFD$FDis,temFD$RaoQ)
colnames(important.indices.T) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
pairs(important.indices.T, pch=20)
# FRis is positively correlated with number of species as expected
# FDis and Rao very similar, actually same index with only a squating difference

data.frame(important.indices.T)%>%
  rownames_to_column(var = "plot")%>%
  gather(FDind, value, NumbSpecies:Rao)%>%
  merge(plot_x_env_T2)%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv"), by = c("plot", "elevation"))%>%
  mutate(FDind = factor(FDind))%>%
  ggplot(aes(y=value, x= elevation, fill = site), color= "black")+ # GDD  FDD Snw  elevation
  geom_point(shape = 21, size =3)+
  geom_smooth (aes(y=value, x= elevation),method = "lm", se=FALSE, color = "black", inherit.aes=F)+
  facet_wrap(~FDind, nrow = 3, scales = "free")+
  labs(title="FD indices in Temperate system")+
  theme_classic(base_size = 12)+
  theme(strip.text = element_text(size =12),
        legend.position = "bottom")

# To test te CWM for microclimatic gradients we could use simple linear model or 
# REML, restricted maximum likelihood model
summary(lm(temFD$FRic~elevation+ FDD + GDD + Snw, data=plot_x_env_T)) # nbsp   RaoQ

# compute the dissimilarity for each trait separately and then average the dissimilarity across traits
# not sure if necessaty with our type of data
head(sp_x_trait_T)
all.dist.T <- (gowdis(sp_x_trait_T["seed_mass"]) + gowdis(sp_x_trait_T["seed_production"]) + 
                 gowdis(sp_x_trait_T["plant_height"]) + gowdis(sp_x_trait_T["floral_height"]) +
                 gowdis(sp_x_trait_T["odds_A_control"]) +gowdis(sp_x_trait_T["odds_B_dark"]) +
                 gowdis(sp_x_trait_T["odds_C_WP"]) +gowdis(sp_x_trait_T["odds_D_constant"]) +
                 gowdis(sp_x_trait_T["odds_E_cold"])) / 9
FD.alltraits.T <- dbFD(all.dist.T, log(plot_x_sp_T+1), message = F, calc.CWM =F, stand.FRic = T)
summary(lm(FD.alltraits.T$FDis ~elevation+ FDD + GDD + Snw, data=plot_x_env_T))

## 5.4.1 cluster dendograma ####
tree.traits.T <- hclust(all.dist.T, "average")
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

## 5.5 Calculation of MDP (main pairwise dissimilarity) and RAo with picante package ####
mpd.T <- mpd(plot_x_sp_T, as.matrix(all.dist.T))
mntd.T <- mntd(plot_x_sp_T, as.matrix(all.dist.T))
number.species.T <- temFD$nbsp

par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
plot(number.species.T, mpd.T, pch = 20, ylab = "MPD",
     xlab = "Number of species")
plot(number.species.T, mntd.T, pch = 20, ylab = "MNTD",
     xlab = "Number of species")
# We observe a general lack of correlation between MPD and the number of species
# as expected, we observe a negative correlation between MNTD and the number of species
## Species abundances can be considered with the argument abundance.weighted = TRUE
# we could logtransforme abundance data to damp a bit the effect of overabundant species

mpd.T.ab <- mpd(log(plot_x_sp_T+1), as.matrix(all.dist.T), abundance.weighted = T)
mntd.T.ab <- mntd(log(plot_x_sp_T+1), as.matrix(all.dist.T), abundance.weighted = T)

plot(FD.alltraits.T$RaoQ, mpd.T.ab, pch=20, ylab= "MPd abundance")

# Use of melodic function to compute both weighted and unweighted forms of MPD and Rao
# activate function from melodic.R script in src folder
melodic.T <- melodic(log(plot_x_sp_T+1), as.matrix(all.dist.T))
str(melodic.T)

# MPD = Rao/Simpson
par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
plot(melodic.T$abundance$mpd, melodic.T$abundance$rao, 
     xlab="MPD", ylab="Rao", pch=20, main= "Temperate community (species rich)")
abline(0, 1)
par(mfrow = c(1, 1))

# Use of Rao function in scr folder
Rao.T <- Rao(log(sp_x_plot_T+1), all.dist.T, dphyl = NULL, weight = F, Jost=T, structure = NULL)

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