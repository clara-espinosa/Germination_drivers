 library (TPD); library (tidyverse); library(ggplot2); library(gridExtra)
# TPD package exploration
# following https://cran.r-project.org/web/packages/TPD/vignettes/TPD.html

#### TPD Species #### 
#Suposse that we want to calculate the TPDs of each species considering 
#Sepal.Length and Sepal.Width. Estimating these TPDs’s is very simple:
traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
sp_iris <- iris$Species
TPDs_iris <- TPDs(species = sp_iris, traits_iris)

#Object of class TPDsp contains:
#the grid of cells in which the function was evaluated, contained in TPDs_iris$data$evaluation_grid, 
#the probability associated to each of those cells, contained in a list TPDs_iris$TPDs 
#with one element for each species:
head(TPDs_iris$data$evaluation_grid)
nrow(TPDs_iris$data$evaluation_grid) 
names(TPDs_iris$TPDs)
head(sort(TPDs_iris$TPDs$setosa, decreasing = TRUE))
sum(TPDs_iris$TPDs$setosa)

#TPDs_iris contains the TPDs of the three species (in TPDs$TPDs), which are 
#vectors containing the probability associated to each point of the evaluation grid
summary(TPDs_iris$TPDs)
sapply(TPDs_iris$TPDs, sum)

#The function plotTPD can be used to visualize the resulting TPDs’s
plotTPD(TPD = TPDs_iris, nRowCol = c(1,3)) 
# Be patient, plots can take some time.
# If you think it takes too long, consider reducing the number of plots using 'whichPlot'
# NOTE: Be careful when making 2 dimensional plots. If you try to plot too many species or 
# populations at the same time, it can take quite a long time. 
#We recommend to use the argument whichPlot in PlotTPD to control for this

# The parameter alpha determines the proportion of the probability density function of each species 
# or population that will be included in the resulting TPD (Blonder et al. 20114; Swanson et al. 2015).
# This means that setting alpha = 1 will include the whole function, whereas alpha = 0.5 will include 
# only a 50%. It is important to remarks that in the cases in which alpha is set to a value smaller 
# than 1, the TPD function rescales the probabilities so that they add up to 1. 
# By default, alpha is set to 0.95; thus TPDs is not too sensitive to outliers.

# Sometimes there is trait information at the popolation level. Imagine that we have measured traits 
# on individuals of each species in two different environments. The argument samples of TPDs allows 
# us to calculate a separate TPDs for each population of each species
pops_iris <- rep(c(rep("Site.1", 25), rep("Site.2", 25)), 3)
TPDs_pops_iris <- TPDs (species = sp_iris, traits = traits_iris, samples = pops_iris)
#Calculating densities for Multiple populations_Multiple species
plotTPD(TPD = TPDs_pops_iris, nRowCol = c(3,2)) 

# The TPDs function requires a high quantity of information to work: one measurement per individual 
# and trait cosidered, and all the traits must be measured in all individuals in order to place 
# them in the multidimensional space. In addition, several individuals are required to estimate 
# reliable TPDs’s

# TPDsMean performs the same task as TPDs, but it estimates the probabilities by calculating a 
# multivariate normal distribution for each species or population, rather than using kernel density 
# estimations. Hence, the only information required is the mean and standard deviation of each 
# trait considered for each species or population

names_iris <- unique(iris$Species)
mt1 <- tapply(iris[, "Sepal.Length"], iris$Species, mean)
mt2 <- tapply(iris[, "Sepal.Width"], iris$Species, mean)
means_iris <- matrix(c(mt1, mt2), ncol=2)
st1 <- tapply(iris[, "Sepal.Length"], iris$Species, sd)
st2 <- tapply(iris[, "Sepal.Width"], iris$Species, sd)
sds_iris <- matrix(c(st1, st2), ncol=2)
TPDsMean_iris<- TPDsMean(species = names_iris, means = means_iris, sds = sds_iris)
#alculating densities for One population_Multiple species
plotTPD(TPD = TPDsMean_iris, nRowCol = c(1,3))

#### TPD Communities #### 
# TPDc combines the information contained in the TPDs of the species in the community with 
# information about their abundances, resulting in the Trait Probabiity Density of the community 
# TPDc is the result of summing the probability of each of the species present in the community,
# weighted by their relative abundances

#Let us imagine three different communities containing different abundances of each species:
abundances_comm_iris <- matrix(c(c(0.5, 0.4, 0), #I. virginica absent
                                 c(0.0, 0.9,  0.1 ), #I. versic. dominates; setosa absent
                                 c(0.0, 0.1,  0.9 )), #I. virg. dominates; setosa absent
                               ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
                                                                       unique(iris$Species))) 

# To calculate the TPDc’s of a group of communities, we first have to calculate the TPDs of all the 
# species present in each community. (This step is important: if you do not have the TPDs of all 
# the species or populations, you cannot estimate the TPDc of the community).
TPDs_iris <- TPDs(species = sp_iris, traits_iris)

TPDc_iris <- TPDc(TPDs = TPDs_iris, sampUnit = abundances_comm_iris )

#Wwe are interested in the element TPDc_iris$TPDc$TPDc, which is a list with an element for each 
# community (length(TPDc_iris$TPDc$TPDc) = 3), each one containing the probability associated to 
# each cell of the grid in each community. 
# In fact, the TPDc of a community reflects the probability of observing a given trait value when 
# extracting a random individual from the community

# plotTPD can be used to provide a graphical representation of TPDc’s of 1 or 2 dimensions.
  plotTPD(TPD = TPDc_iris, nRowCol = c(1,3))

## 1. Species are not of primary concern
### 1.a. Indices for single units #### 
# The function REND calculates these three indices using the TPDs’s and TPDc’s of populations/species
# and communities, respectively. The function works with one or more dimensions, thus expanding 
# the framework presented in Mason et al. (2005) to multiple traits without modifying its original 
# rationale, based on probabilities. Let us start with species with one and two dimensions:

TPDs_iris_d2 <- TPDs(species = sp_iris, traits_iris)
RED_d2 <- REND(TPDs = TPDs_iris_d2)
# RED_d2 contain the FRic, FEve and FDiv values of the three species considering 2 traits simultaneously. 
# We can plot the respective TPDs’s, along with their FRich values
plotTPD(TPD = TPDs_iris_d2, nRowCol = c(1, 3), leg.text = paste0(names(RED_d2$species$FRichnes),
                                                                 "; FRic=", round(RED_d2$species$FRichnes, 3)))
#The first thing to underscore is that FRic is given in the original trait units. FRic is expressed in cm and 
# in the two dimensional case (Sepal.Length and Sepal.Width) it is expressed in cm2
# The plots also reveal that I. setosa occupies a portion of the functional space that is distinct 
# from the ones occupied by the other two species

# We are going to take advantage of this observation to show the calculation of FDiv, and 
# its application to communities
abundances_comm_iris_divergence <- matrix(c(c(0.45, 0.1, 0.45), # Different species are more abundant
                                            c(0.05, 0.9,  0.05 ),  # Different species are less abundant
                                            c(0.9,  0.05, 0.05 )), # A species with extreme trait values is very abundant
                                          ncol = 3, byrow = TRUE, dimnames = list(paste0("Comm.",1:3),
                                                                                  unique(iris$Species))) 
TPDc_iris_divergence_d2 <- TPDc(TPDs = TPDs_iris_d2, sampUnit = abundances_comm_iris_divergence)
RED_comm_d2 <- REND(TPDc = TPDc_iris_divergence_d2)

plotTPD(TPD = TPDc_iris_divergence_d2, nRowCol = c(1, 3), 
        leg.text = paste0(names(RED_comm_d2$communities$FDivergence),
                          "; FDiv=", round(RED_comm_d2$communities$FDivergence, 3)))
# Comm.3, with a great proportion of their TPDc being away from the center, displays the highest 
# FDiv value, followed by Comm.1
### 1.b. Indices for multiple units ####
# Overlap-based functional dissimilarity, and its decomposition (dissim)
plotTPD(TPDc_iris, nRowCol = c(1,3))
dissim_iris_comm <- dissim(TPDc_iris)
dissim_iris_comm$communities$dissimilarity
dissim_iris_comm$communities$P_shared
dissim_iris_comm$communities$P_non_shared

