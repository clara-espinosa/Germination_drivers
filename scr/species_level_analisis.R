library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(vegan); library(ade4); library(MASS)
library(ecodist);library(maptools);library(rpart);library(splines)
library(gam);library(pgirmess);library(utils);library(combinat)
library(cluster);library(fpc);library(clusterSim);library(lmtest)
library(Hmisc);library(gplots);library(NbClust);library(rpart)
library(rpart.plot);library(dismo);library(multcomp);library(gbm)
library(raster);library(tidyverse);
library(glmmTMB); library (DHARMa)

# following chapter 4 multivariate species level responses by Lars ###
# Separate from the beggining between communities
############################################# MEDITERRANEAN ###########################################################
# 4.3 Trait free CCA ####
cca1M <- cca(sp_x_plot_M ~ elevation+GDD+FDD, data=plot_x_env_M) # +bio1+bio2+bio7+Snw tried but highly correlated
plot(cca1M, display = c("species", "bp"))
anova(cca1M, by= "terms") # elevation and GDD significantly explain the ordination of the species
RsquareAdj(cca1M) # roughly 6% of variation in the species data is explained by these variables
head(vegan::scores(cca1M, display = "species"))
# possibility to separate CCA for each env. variable and then use the first CCA axes to obtain 
# the value (score) for the particular gradient

# 4.4 Double canonical correspondance analysis (dbrda function not working!!) ####
ca1M <- dudi.coa(sp_x_plot_M, scannf = F)
dCCA1M <- dbrda(ca1M, plot_x_env_M, sp_x_trait_M, scannf = FALSE)
# 4.5 Functional response groups ####
dist_CCAM <- dist(vegan::scores(cca1M, display = "species"))
clust_CCAM <- hclust(dist_CCAM, method = "ward.D2")
plot(clust_CCAM, cex = 0.6)

#The function NbClust from the package of the same name estimates the optimal number of groups on the basis of maximizing 
#the dissimilarity between groups and minimizing the dissimilarity within groups.
groups_CCAM <- NbClust(diss = dist_CCAM, distance = NULL, min.nc = 2, max.nc = 6,
                       method = "ward.D2", index = "silhouette")
groups_CCAM # optimal 3 clusters

# we are asking whether there are traits that explain to which cluster the species belong, keeping in
# mind that this clustering is based on our initial ordination that arranges species along the
# considered environmental gradients. For this test we use here a simple regression analysis
# test every trait against the cluster affiliation
summary(lm(as.matrix(sp_x_trait_M) ~ groups_CCAM$Best.partition))
# No trait significantly associated to each of the groups

# 4.6 RDA and regression trees (more appropiate with short gradients and sp abundance-gradient linear relationship) ####
rda1M <- rda(sp_x_plot_M ~ ., data = plot_x_env_M, scale = TRUE)
plot(rda1M, display = c("bp", "sp"))
RsquareAdj(rda1M)
#how much of the variation the 3 constraining factors explain.
cumsum(rda1M$CCA$eig) / sum(rda1M$CCA$eig)

# Analyzing the relationship between the species responses (expressed as the species scores in
#the RDA ordination space) and their traits
lm_rda1M <- lm(vegan::scores(rda1M, choices = 1, display = "species") ~ ., data = sp_x_trait_M)
summary(lm_rda1M)
# significance of odds_D_constant
# marginal significance seed mass, floral height, D_constant_proportion
# the above traits of our species do no explain any significant proportion of the sp responses
### STOPED FOLLOWING CHAPTER 4 HERE, NEXT WOULD BE REGRESSION TREES

################################# TEMPERATE ########################################################################
# 4.3 Trait free CCA ####
cca1T <- cca(sp_x_plot_T ~ elevation+GDD+FDD+bio1+bio2+bio7+Snw, data=plot_x_env_T) #  tried but highly correlated
plot(cca1T, display = c("species", "bp"))
anova(cca1T, by= "terms") # elevation, GDD, FDD, bio7 and snow significantly explain the ordination of the species
RsquareAdj(cca1T) # roughly 15% of variation in the species data is explained by these variables
head(vegan::scores(cca1T, display = "species"))
# separate CCA for each env. variable and then use the first CCA axes to obtain the value for the particular gradient
# still needs to account for the effects of the other variables as conditional variables (covariates)

# 4.4 Double canonical correspondance analysis (dbrda function not working!!) ####
ca1T<- dudi.coa(sp_x_plot_T, scannf = F)
dCCA1T <- dbrda(ca1T, plot_x_env_T, sp_x_trait_T, scannf = FALSE)

# 4.5 Functional response groups ####
dist_CCAT <- dist(vegan::scores(cca1T, display = "species"))
clust_CCAT <- hclust(dist_CCAT, method = "ward.D2")
plot(clust_CCAT, cex = 0.6)

#The function NbClust from the package of the same name estimates the optimal number of groups on the basis of maximizing 
#the dissimilarity between groups and minimizing the dissimilarity within groups.
groups_CCAT <- NbClust(diss = dist_CCAT, distance = NULL, min.nc = 2, max.nc = 6,
                       method = "ward.D2", index = "silhouette")
groups_CCAT# optimal 2 clusters

# we are asking whether there are traits that explain to which cluster the species belong, keeping in
# mind that this clustering is based on our initial ordination that arranges species along the
# considered environmental gradients. For this test we use here a simple regression analysis
# test every trait against the cluster affiliation
summary(lm(as.matrix(sp_x_trait_T) ~ groups_CCAT$Best.partition))
# No seed trait significantly associated to each of the groups

# 4.6 RDA and regression trees (more appropiate with short gradients and sp abundance-gradient linear relationship) ####
rda1T <- rda(sp_x_plot_T ~ ., data = plot_x_env_T, scale = TRUE)
plot(rda1T, display = c("bp", "sp"))
RsquareAdj(rda1T)
#how much of the variation the 3 constraining factors explain.
cumsum(rda1T$CCA$eig) / sum(rda1T$CCA$eig)

# Analyzing the relationship between the species responses (expressed as the species scores in
#the RDA ordination space) and their traits
lm_rda1T <- lm(vegan::scores(rda1T, choices = 1, display = "species") ~ ., data = sp_x_trait_T)
summary(lm_rda1T)
# significance of seed mass, plant height, odds_A_control