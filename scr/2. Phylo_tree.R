library(tidyverse);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)
library (tidyverse); library(scales)

### Phylo tree both communities #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =",") %>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
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
#devtools::install_github("jinyizju/V.PhyloMaker2")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree.tree")
# both communities tree without species wth 0 germ necessary??####
read.csv("data/species.csv", sep =",") %>%
  mutate(species= str_replace(species, "Minuartia sp", "Minuartia arctica"))%>%
  select (species, family) %>%
  unique %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianella campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea") %>%
  filter(!species=="Phalacrocarpum oppositifolium")%>% # 0 germ across all treatments
  filter(!species == "Teesdalia conferta")%>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1

phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/treegerm.tree")
# we already have the phylotree for both communities made in phylo_tree script
# with all species and also without the species that had 0 germination
