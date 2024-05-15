library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library(ape);library(phytools);library(Rphylopars);library(V.PhyloMaker);library(doBy)

## TRAITS AND PHYLOGENY (Chapter 8 Lars book) ###
par(mfrow = c(1, 1))
## Phylogenetic trees ####
# we already have the phylotree for both communities made in phylo_tree script
# with all species and also without the species that had 0 germination
# make phylo tree with phylomaker for each of our communitites
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
#mediterranean
read.csv("data/species.csv") %>%
  filter(community == "Mediterranean")%>%
  select (species, family) %>%
  unique %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1
unique(ranks1$species)
#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/MED_tree.tree")

#temperate
read.csv("data/species.csv") %>%
  filter(community == "Temperate")%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  select (species, family) %>%
  unique %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1
unique(ranks1$species)
#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/TEM_tree.tree")


###########################################  MEDITERRANEAN  ##################################################
## 8.4. Trait phylogenetic signal ####
# dataframe with trait x species #
str(sp_x_trait_M) 
rownames(sp_x_trait_M) # 24 species as rownames
# dataframe with phylogenetic data
MED_tree <- read.tree("results/MED_tree.tree")
plot(MED_tree)
MED_tree$tip.label <- gsub("_", " ", MED_tree$tip.label) # 24 tips labels

# match data function from Lars book in match_data.R script in src folder
match_traits_phylo_M <- match_data(traits=na.omit(sp_x_trait_M), tree = MED_tree) # here seems to work properly
#$tree says "Phylogenetic tree with 24 tips and 22 internal nodes2
str(match_traits_phylo_M) # $tree : List of 5 only matches 5 first species???
head(match_traits_phylo_M$tree$tip.label) # complete list with tip labels
plot(match_traits_phylo_M$tree) ## tree looks okei when plotted

#extract each single trait variables
#seed mass
sm <- as.matrix(match_traits_phylo_M$traits[ ,1])
phytools::phylosig (match_traits_phylo_M$tree, sm, method = "K") # 0.45
phytools::phylosig (match_traits_phylo_M$tree, sm, method = "lambda") ##0.0437 logL(lambda) = -132
#seed production
sp <- as.matrix(match_traits_phylo_M$traits[ ,2])
phytools::phylosig (match_traits_phylo_M$tree, sp, method = "K") # 0.92
phytools::phylosig (match_traits_phylo_M$tree, sp, method = "lambda") ##1.05logL(lambda) = -148
#plant height
ph <- as.matrix(match_traits_phylo_M$traits[ ,3])
phytools::phylosig (match_traits_phylo_M$tree, ph, method = "K") # 0.33
phytools::phylosig (match_traits_phylo_M$tree, ph, method = "lambda") ##0.00008logL(lambda) = -88
#floral height
fh <- as.matrix(match_traits_phylo_M$traits[ ,4])
phytools::phylosig (match_traits_phylo_M$tree, fh, method = "K") # 0.26
phytools::phylosig (match_traits_phylo_M$tree, fh, method = "lambda") ##0.19 logL(lambda) = -92
# odds_dark
od <- as.matrix(match_traits_phylo_M$traits[ ,5])
phytools::phylosig (match_traits_phylo_M$tree, od, method = "K") # 0.57
phytools::phylosig (match_traits_phylo_M$tree, od, method = "lambda") ##1.01logL(lambda) = -90
# odds_WP
ow <- as.matrix(match_traits_phylo_M$traits[ ,6])
phytools::phylosig (match_traits_phylo_M$tree, ow, method = "K") # 0.42
phytools::phylosig (match_traits_phylo_M$tree, ow, method = "lambda") ##0.00008logL(lambda) = -92
# odds_constant
otc <- as.matrix(match_traits_phylo_M$traits[ ,7])
phytools::phylosig (match_traits_phylo_M$tree, otc, method = "K") # 0.36
phytools::phylosig (match_traits_phylo_M$tree, otc, method = "lambda") ##0.00008 logL(lambda) = -89
# odds_cold
oc <- as.matrix(match_traits_phylo_M$traits[ ,8])
phytools::phylosig (match_traits_phylo_M$tree, oc, method = "K") # 0.27
phytools::phylosig (match_traits_phylo_M$tree, oc, method = "lambda") ##0.09 logL(lambda) = -80


## 8.5 Traitgrams  ####
# test with only 2 plots out of 73 of the mediterranean community
plot1_M <- plot_x_sp_M[1,] [,plot_x_sp_M[1,]!=0] # second half of the command !=0 not working
# like this seems to work
plot1_M <- plot_x_sp_M[1,(plot_x_sp_M[1,]!=0)] 
plot1_M <- plot_x_sp_M[,(plot_x_sp_M[1,]!=0)] # like this match_data function works
plot2_M <- plot_x_sp_M[6,(plot_x_sp_M[6,]!=0)] 
plot2_M <- plot_x_sp_M[,(plot_x_sp_M[6,]!=0)]# like this match_data function works
# match_data function to join all data frames (not working)
# when only considering 1 row (plot) drop all tips of the tree

match_plot1_M <- match_data(com= as.matrix(plot1_M), traits = na.omit(sp_x_trait_M), tree=MED_tree)
match_plot2_M <- match_data(com= as.matrix(plot2_M), traits = na.omit(sp_x_trait_M), tree=MED_tree)

par(mfrow = c(1, 2))
traitgram(as.matrix(match_plot1_M$traits["seed_mass"]), multi2di(match_plot1_M$tree),
          show.names = TRUE, cex = 0.7)
# PROBLEM "trait and phy names do not match"
traitgram(as.matrix(match_plot2_M$traits["seed_mass"]), multi2di(match_plot2_M$tree),
          show.names = TRUE)

## 8.6 Phylogenetic comparative methods (at species level) only done for few traits ####
tree_pic_M <- multi2di(match_traits_phylo_M$tree)
seedmass_M <- match_traits_phylo_M$traits[,"seed_mass"]
seedprod_M <- match_traits_phylo_M$traits[,"seed_production"]
Pheight_M <- match_traits_phylo_M$traits[,"plant_height"]
Fheight_M <- match_traits_phylo_M$traits[,"floral_height"]
odds_B_dark_M <- match_traits_phylo_M$traits[,"odds_B_dark"]
odds_C_WP_M <- match_traits_phylo_M$traits[,"odds_C_WP"]
odds_D_constant_M <- match_traits_phylo_M$traits[,"odds_D_constant"]
odds_E_cold_M <- match_traits_phylo_M$traits[,"odds_E_cold"]
# We use the "pic" function from ape package
pic_seedmass_M<- pic(seedmass_M, tree_pic_M)
pic_seedprod_M<- pic(seedprod_M, tree_pic_M)
pic_Pheight_M<- pic(Pheight_M, tree_pic_M)
pic_Fheight_M<- pic(Fheight_M, tree_pic_M)
pic_odds_B_dark_M<- pic(odds_B_dark_M, tree_pic_M)
pic_odds_C_WP_M<- pic(odds_C_WP_M, tree_pic_M)
pic_odds_D_constant_M<- pic(odds_D_constant_M, tree_pic_M)
pic_odds_E_cold_M<- pic(odds_E_cold_M, tree_pic_M)

# the relationship between PICS can be analysed using a linear model, fitting intercept through the origin (-1)
lm_pics <- lm(pic_seedmass_M~pic_seedprod_M-1 )
summary(lm_pics)
lm_ord <- lm (seedmass_M~seedprod_M)
summary(lm_ord)
# without taking into account phylogeny significant, with pics NS
par(mfrow=c(1,2))
plot(pic_seedmass_M, pic_seedprod_M)
abline(lm_pics, col="orange")
plot(seedmass_M,seedprod_M)
abline(lm_ord, col="orange")
par(mfrow=c(1,1))
## 8.7 Imputing trait data with the help of phylogeny ####
# For now all species with traits of interest
## 8.8 Phylogenetic diversity  (melodic function giving ERROR ####
match_PD_M <- match_data (com = plot_x_sp_M, tree= MED_tree)
phylo_dist_M <- cophenetic(match_PD_M$tree)
phylo_dist_M

# We show the main pairwise distance index function mpd but with the version from the melodic function
# this function expects the distance matrix to be an actual amtrix object which is already the case for 
#phylogenetic distances produced by cophenetic function
melody_M <- melodic(match_PD_M$com, phylo_dist_M, type = "abundance")
#Error in if (class(samp) != "matrix") { : the condition has length > 1

## How to combine trait and phylogenetic diversity ####
dc_M <- decouple(traits = na.omit(sp_x_trait_M[,"seed_mass", drop = FALSE]), tree=MED_tree)
lapply(dc_M, class)
###########################################  TEMPERATE ##################################################
## 8.4. Trait phylogenetic signal ####
# dataframe with trait x species #
str(sp_x_trait_T) 
rownames(sp_x_trait_T) # 24 species as rownames
# dataframe with phylogenetic data
TEM_tree <- read.tree("results/TEM_tree.tree")
plot(TEM_tree)
TEM_tree$tip.label <- gsub("_", " ", TEM_tree$tip.label) # 24 tips labels

# match data function from Lars book in match_data.R script in src folder
match_traits_phylo_T <- match_data(traits=na.omit(sp_x_trait_T), tree = TEM_tree) # here seems to work properly
#$tree says "Phylogenetic tree with 24 tips and 22 internal nodes2
str(match_traits_phylo_T) # $tree : List of 5 only matches 5 first species???
head(match_traits_phylo_T$tree$tip.label) # complete list with tip labels
plot(match_traits_phylo_T$tree) ## tree looks okei when plotted

#extract each single trait variables
#seed mass
sm <- as.matrix(match_traits_phylo_T$traits[ ,1])
phytools::phylosig (match_traits_phylo_T$tree, sm, method = "K") # 1.1
phytools::phylosig (match_traits_phylo_T$tree, sm, method = "lambda") ##1.02 logL(lambda) = -158
#seed production
sp <- as.matrix(match_traits_phylo_T$traits[ ,2])
phytools::phylosig (match_traits_phylo_T$tree, sp, method = "K") # 0.34
phytools::phylosig (match_traits_phylo_T$tree, sp, method = "lambda") ##0.00007logL(lambda) = -185
#plant height
ph <- as.matrix(match_traits_phylo_T$traits[ ,3])
phytools::phylosig (match_traits_phylo_T$tree, ph, method = "K") # 0.28
phytools::phylosig (match_traits_phylo_T$tree, ph, method = "lambda") ##0.26logL(lambda) = -78
#floral height
fh <- as.matrix(match_traits_phylo_T$traits[ ,4])
phytools::phylosig (match_traits_phylo_T$tree, fh, method = "K") # 0.23
phytools::phylosig (match_traits_phylo_T$tree, fh, method = "lambda") ##0.27 logL(lambda) = -98
# odds_dark
od <- as.matrix(match_traits_phylo_T$traits[ ,5])
phytools::phylosig (match_traits_phylo_T$tree, od, method = "K") # 0.31
phytools::phylosig (match_traits_phylo_T$tree, od, method = "lambda") ##0.16 logL(lambda) = -135
# odds_WP
ow <- as.matrix(match_traits_phylo_T$traits[ ,6])
phytools::phylosig (match_traits_phylo_T$tree, ow, method = "K") # 0.28
phytools::phylosig (match_traits_phylo_T$tree, ow, method = "lambda") ##0.00007logL(lambda) = -132
# odds_constant
otc <- as.matrix(match_traits_phylo_T$traits[ ,7])
phytools::phylosig (match_traits_phylo_T$tree, otc, method = "K") # 0.32
phytools::phylosig (match_traits_phylo_T$tree, otc, method = "lambda") ##0.01 logL(lambda) = -129
# odds_cold
oc <- as.matrix(match_traits_phylo_T$traits[ ,8])
phytools::phylosig (match_traits_phylo_T$tree, oc, method = "K") # 0.38
phytools::phylosig (match_traits_phylo_T$tree, oc, method = "lambda") ##0.00007 logL(lambda) = -126


## 8.5 Traitgrams  ####
# test with only 2 plots out of 73 of the mediterranean community
plot1_T <- plot_x_sp_T[1,] [,plot_x_sp_T[1,]!=0] # second half of the command !=0 not working
# like this seems to work
plot1_T <- plot_x_sp_T[1,(plot_x_sp_T[1,]!=0)] 
plot1_T <- plot_x_sp_T[,(plot_x_sp_T[1,]!=0)] # like this match_data function works
plot2_T <- plot_x_sp_T[6,(plot_x_sp_T[6,]!=0)] 
plot2_T <- plot_x_sp_T[,(plot_x_sp_T[6,]!=0)]# like this match_data function works
# match_data function to join all data frames (not working)
# when only considering 1 row (plot) drop all tips of the tree

match_plot1_T <- match_data(com= as.matrix(plot1_T), traits = na.omit(sp_x_trait_T), tree=TEM_tree)
match_plot2_T <- match_data(com= as.matrix(plot2_T), traits = na.omit(sp_x_trait_T), tree=TEM_tree)

par(mfrow = c(1, 2))
traitgram(as.matrix(match_plot1_T$traits["seed_mass"]), multi2di(match_plot1_T$tree),
          show.names = TRUE, cex = 0.7)
# PROBLEM "trait and phy names do not match"
traitgram(as.matrix(match_plot2_T$traits["seed_mass"]), multi2di(match_plot2_T$tree),
          show.names = TRUE)

## 8.6 Phylogenetic comparative methods (at species level) only done for few traits ####
tree_pic_T <- multi2di(match_traits_phylo_T$tree)
seedmass_T <- match_traits_phylo_T$traits[,"seed_mass"]
seedprod_T <- match_traits_phylo_T$traits[,"seed_production"]
Pheight_T <- match_traits_phylo_T$traits[,"plant_height"]
Fheight_T <- match_traits_phylo_T$traits[,"floral_height"]
odds_B_dark_T <- match_traits_phylo_T$traits[,"odds_B_dark"]
odds_C_WP_T <- match_traits_phylo_T$traits[,"odds_C_WP"]
odds_D_constant_T <- match_traits_phylo_T$traits[,"odds_D_constant"]
odds_E_cold_T <- match_traits_phylo_T$traits[,"odds_E_cold"]
# We use the "pic" function from ape package
pic_seedmass_T<- pic(seedmass_T, tree_pic_T)
pic_seedprod_T<- pic(seedprod_T, tree_pic_T)
pic_Pheight_T<- pic(Pheight_T, tree_pic_T)
pic_Fheight_T<- pic(Fheight_T, tree_pic_T)
pic_odds_B_dark_T<- pic(odds_B_dark_T, tree_pic_T)
pic_odds_C_WP_T<- pic(odds_C_WP_T, tree_pic_T)
pic_odds_D_constant_T<- pic(odds_D_constant_T, tree_pic_T)
pic_odds_E_cold_T<- pic(odds_E_cold_T, tree_pic_T)

# the relationship between PICS can be analysed using a linear model, fitting intercept through the origin (-1)
lm_pics <- lm(pic_seedmass_T~pic_seedprod_T-1 )
summary(lm_pics)
lm_ord <- lm (seedmass_T~seedprod_T)
summary(lm_ord)
# without taking into account phylogeny significant, with pics NS
par(mfrow=c(1,2))
plot(pic_seedmass_T, pic_seedprod_T)
abline(lm_pics, col="orange")
plot(seedmass_T,seedprod_T)
abline(lm_ord, col="orange")
par(mfrow=c(1,1))
## 8.7 Imputing trait data with the help of phylogeny ####
# For now all species with traits of interest
## 8.8 Phylogenetic diversity  (melodic function giving ERROR ####
match_PD_T <- match_data (com = plot_x_sp_T, tree= TEM_tree)
phylo_dist_T <- cophenetic(match_PD_T$tree)
phylo_dist_T

# We show the main pairwise distance index function mpd but with the version from the melodic function
# this function expects the distance matrix to be an actual amtrix object which is already the case for 
#phylogenetic distances produced by cophenetic function
melody_T <- melodic(match_PD_T$com, phylo_dist_T, type = "abundance")
#Error in if (class(samp) != "matrix") { : the condition has length > 1

## How to combine trait and phylogenetic diversity ####
dc_T <- decouple(traits = na.omit(sp_x_trait_T[,"seed_mass", drop = FALSE]), tree=TEM_tree)
lapply(dc_M, class)
