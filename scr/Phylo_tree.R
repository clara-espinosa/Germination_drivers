library(tidyverse);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)
library (tidyverse); library(scales)

### Phylo tree both communitites #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =",") %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album CF", "Sedum album")) %>%
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

#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree.tree")
