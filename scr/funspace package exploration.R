library(tidyverse); library (funspace)

################################# MEDITERRANEAN ##############################
# complete trait data
germ_odds_M%>% # log odds ratios
  merge(seed_mass_M)%>% # log transformed
  merge(plant_height_M)%>% # log transformed
  merge(leaves_traits_M )%>% #only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height,leaf_area, LDMC, SLA,
                odds_B_dark, odds_C_WP, odds_D_constant)%>% 
  column_to_rownames(var="species")-> M.data.funspace
##### Module 1: Building and exploring the functional trait space ####
# testing the number of dimensions defining the functional trait space
funspaceDim(M.data.funspace)

# Run PCA 
M.pca.funspace <- princomp(M.data.funspace, cor=T)

# building the functional trait space (using the first two Pcs)
M.trait_space_global <- funspace(x = M.pca.funspace, PCs = c(1,2), n_divisions = 300)
summary(M.trait_space_global)
# first axis mainly explaining leaf area, seed mass but also dark and constant T
# second axis mainly explaining water stress, SLA and LDMC
# Overall explained: quality of the representation of the functional space
    # leaf_area 74.9, water stress 71.7, LDMC 66.7, seed mass 61.7
# Functional diversity indicators:
    # functional richness 99.52; functional divergence = 0.45

##### Module 2. Mapping a functional trait space ####
# Response variable we want to map
##### Module 3. Plotting a functional trait space ####
plot(M.trait_space_global, quant.plot=T, 
     pnt = T, pnt.cex = 1,pnt.col = rgb(0.2, 0.8, 0.1, alpha = 0.2),
     arrows = TRUE, ,arrows.length = 0.9, 
     title(main = "PCA Mediterranean community"))

################################ TEMPERATE ###########################
germ_odds_T%>% # log odd ratios
  merge(seed_mass_T)%>% #log seed mass
  merge(plant_height_T)%>% # log plant height
  merge(leaves_traits_T)%>% # only leaf area log transformed
  dplyr::select(species, seed_mass,plant_height, 
                leaf_area, LDMC, SLA, 
                odds_B_dark, odds_C_WP, odds_D_constant)%>%
  column_to_rownames(var="species")->T.data.funspace
##### Module 1: Building and exploring the functional trait space ####
# testing the number of dimensions defining the functional trait space
funspaceDim(T.data.funspace)

# Run PCA 
T.pca.funspace <- princomp(T.data.funspace, cor=T)

# building the functional trait space (using the first two Pcs)
T.trait_space_global <- funspace(x = T.pca.funspace, PCs = c(1,2), n_divisions = 300)
summary(T.trait_space_global)
# first axis mainly explaining leaf area and plant height
# second axis mainly explaining constant T, water stress and darkness
# Overall explained: quality of the representation of the functional space
# constant Temp 85.57, leaf_area 75.09, plant height 64.37 
# Functional diversity indicators:
# functional richness 84.93; functional divergence = 0.42

##### Module 2. Mapping a functional trait space ####
# Response variable we want to map
##### Module 3. Plotting a functional trait space ####
plot(T.trait_space_global, quant.plot=T, 
     pnt = T, pnt.cex = 1,pnt.col = rgb(0.2, 0.8, 0.1, alpha = 0.2),
     arrows = TRUE,arrows.length = 0.9, 
     title(main = "PCA Temperate community"))
