library(tidyverse);library(readxl);library(ggplot2);library(dplyr)
library(viridis);library(rstatix);library(lubridate);library(vegan)
library(rstatix);library(geomtextpath);library(psych)
library (ggrepel);library (ggpubr) ;library(binom)
library(RColorBrewer)

# treatments colors ####
brewer.pal(n = 5, name = "Spectral")
"#D7191C" "#FDAE61" "#FFFFBF" "#ABDDA4" "#2B83BA"

######## EXTRA PLOTS playing with VISUALIZATION #####
# germination curves / rate  x community and x treatment (update with error bars) #### 
x11()
germ_df %>%
  mutate ( D0 = 0) %>%
  dplyr::select (species, treatment, viable, D0, D7, D14, D21, D28, D35, D42) %>%
  gather ("scores", "germ", D0:D42) %>%
  filter (!(treatment == "alternate_dark")) %>%
  mutate (germ = replace_na(germ, 0)) %>%
  mutate (scores = factor(scores, levels = c("D0","D7", "D14", "D21", "D28","D35", "D42"))) %>%
  arrange(species, treatment, scores)%>% 
  merge (species) %>%
  filter (!species == "Euphrasia salisburgensis")%>%
  filter (!species == "Gentiana verna")%>%
  filter (!species == "Gentianela campestris")%>%
  filter (!species == "Kobresia myosuroides")%>%
  filter (!species == "Salix breviserrata")%>%
  filter (!species == "Sedum album")%>%
  filter (!species == "Sedum atratum")%>%
  filter (!species == "Solidago virgaurea")%>%
  filter (!species == "Teesdalia conferta")%>%
  filter (!species == "Phalacrocarpum oppositifolium")%>%
  group_by (community, treatment, scores)%>%
  summarise (germinated = sum(germ), viable = sum(viable))%>%
  mutate (cum_germ = cumsum(germinated))  %>% #
  mutate(germpro = (cum_germ/viable), # 
         germpro = round (germpro, digit =2)) %>%
  filter(!treatment == "B_alternate_dark")%>%
  ggplot(aes(x= scores, y= germpro, color= treatment, group = treatment))+
  geom_line(linewidth = 1.5) +
  facet_grid(~community) +
  ylim (0,1)+
  scale_color_manual (name= "Treatment", labels = c("Control", "Water stress", "Constant Temperature" ), 
                      values = c("#2B83BA",  "#FFFFBF","#FDAE61")) +
  #scale_color_viridis_d() +
  labs (title = "Cumulative Germination")+
  xlab ("Germination scores") +
  ylab("Germination proportion") +
  theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         legend.position = "right")

#loop for many graphs   

for (var in unique(cumgerm_sp_treatment$species)) {
  treat_plot = ggplot(cumgerm_sp_treatment[cumgerm_sp_treatment$species==var,], 
                      aes(x= scores, y= germpro, color= treatment, group = treatment))+
    geom_line(linewidth = 2) +
    ylim (0,1)+
    scale_color_manual (name= "Treatment", values = c ( "A_alternate_light" = "#440154FF", "B_alternate_dark" = "#2A788EFF","C_alternate_WP"= "#5ec962" ,"D_constant_light" = "#fde725")) +
    labs (title= var, x= "Germination scores", y ="Germination proportion")+
    theme_classic (base_size = 16) + #theme_minimal for all species for mean treatment
    theme (plot.title = element_text ( size = 32), #hjust = 0.5,
           legend.position = "right")
  theme(plot.title = element_text (size = 30),
        #strip.text = element_text (size = 24, face = "italic"),
        axis.title.y = element_text (size=24), 
        axis.title.x = element_text (size=24), 
        axis.text.x= element_text (size=18, vjust = 0.5),
        legend.position = "none",
        plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm"))
  ggsave(treat_plot, file = paste0("Germination ", var,".png"),
         path = NULL, scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}


# germination curves / rate  x species and x treatment IF USED NEEDS TO BE UPDATED #### 
# cold strat x species ordered by phylogeny ####
  scale_fill_manual(values = c("#F2F9AA"= "Asteraceae", "#F3ec70"="Campanulaceae","#f6e528"="Apiaceae",
                               "#Fad220" = "Lamiaceae","#fcb622" = "Orobanchaceae","darkgoldenrod3" ="Plantaginaceae",
                               "#f68b08" = "Gentianaceae","#ff7125" ="Primulaceae","#c87107" ="Caryophyllaceae",
                               "#8E5005" ="Plumbaginaceae","tomato2"="Rubiaceae",
                               "#08f9dd" ="Fabaceae","#16cacb" ="Salicaceae","#21a8be"="Brassicaceae",
                               "#2d7faf" ="Cistaceae", "#275381"="Crassulaceae","#42346f"="Saxifragaceae",
                               "#45A747" ="Cyperaceae","#396F3E"= "Juncaceae", "#78E208" ="Poaceae","springgreen2"= "Asparagaceae"))
  
  col_germ <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220","#fcb622", 
                "darkgoldenrod3", "#ff7125","#c87107","#8E5005","tomato2", 
                "#08f9dd",  "#21a8be", "#2d7faf", "#275381", "#42346f",
                "#78E208", "#45A747", "#396F3E", "springgreen2")
  x11()
  cold_strat %>%
    filter (!species == "Euphrasia salisburgensis")%>% # specialist
    filter (!species == "Gentiana verna")%>% # specialist
    filter (!species == "Gentianella campestris")%>% # generalist
    filter (!species == "Kobresia myosuroides")%>% #specialist
    filter (!species == "Salix breviserrata")%>% # specialist
    filter (!species == "Sedum album")%>% # generalist
    filter (!species == "Sedum atratum")%>% # generalist
    mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
    merge (species)%>%
    dplyr::select (community, family, species, code,  treatment, finalgerm, viable, mean, upper, lower)%>%
    mutate(family = as.factor(family))%>%
    mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae",
                                "Plantaginaceae","Primulaceae","Caryophyllaceae", "Plumbaginaceae", "Rubiaceae", 
                                "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                                "Crassulaceae", "Saxifragaceae", "Poaceae", 
                                "Cyperaceae", "Juncaceae", "Asparagaceae"))%>%
    mutate(species = as.factor(species))%>%
    mutate(species = fct_relevel(species,"Scilla verna","Luzula caespitosa" , "Carex sempervirens",
                                 "Koeleria vallesiana","Helictochloa marginata", "Festuca glacialis", "Agrostis tileni",  
                                 "Festuca summilusitana","Festuca rubra", "Avenella flexuosa","Neoschischkinia truncatula", 
                                 "Saxifraga paniculata", "Saxifraga oppositifolia" ,"Saxifraga conifera",
                                 "Sempervivum arachnoideum","Sempervivum vicentei","Sedum brevifolium", "Sedum anglicum", 
                                 "Helianthemum canum", "Helianthemum urrielense",
                                 "Teesdalia conferta","Anthyllis vulneraria", "Galium pyrenaicum", 
                                 "Armeria cantabrica", "Armeria duriaei", "Spergula morisonii", "Minuartia CF",
                                 "Minuartia recurva", "Minuartia verna", "Cerastium ramosissimum",
                                 "Arenaria erinacea", "Gypsophila repens", "Dianthus langeanus", 
                                 "Silene ciliata", "Silene acaulis","Androsace villosa",
                                 "Plantago alpina", "Plantago holosteum", "Linaria alpina",
                                 "Veronica nummularia", "Pedicularis pyrenaica", "Thymus praecox", 
                                 "Conopodium majus", "Jasione cavanillesii", "Jasione laevis",
                                 "Phyteuma hemisphaericum", "Jurinea humilis", "Solidago virgaurea", 
                                 "Phalacrocarpum oppositifolium", "Erigeron alpinus"))%>%
    group_by (community, family, species)%>%
    summarize(germpro = mean(mean), 
              upper = mean(upper),
              lower = mean (lower)) %>% 
    ggplot(aes(x= species, y= germpro, fill = family), color = "black") +
    geom_bar (stat = "identity") +
    scale_fill_manual(values = col_germ)+
    geom_errorbar( aes(species, germpro, ymin = lower, ymax = upper), width = 0.5, linewidth = 0.5) +
    coord_flip() + 
    labs (title = "A) Germination during cold stratification", x = "Species", y = "Germination proportion")+
    theme_classic (base_size = 12) +
    theme (plot.title = element_text ( size = 14), #hjust = 0.5,
           legend.position = "none")
  


# comparison between optimals temperatures
  finalgerm%>%
    merge(read.csv("data/species.csv"), by = c("species", "code"))%>%
    filter(!treatment== "E_cold_stratification")%>%
    group_by(treatment, opt_temp)%>%
    summarise(finalgerm=sum(finalgerm),
              viable=sum(viable))%>%
    mutate (binom.confint(finalgerm, viable, methods = "wilson"))%>%
    ggplot(aes(x= treatment, y= mean, fill= opt_temp))+
    geom_bar(color = "black", stat = "identity",position = "dodge", width = 0.9) +# 
    geom_errorbar(aes(x= treatment, y=mean, ymin=lower, ymax=upper), color = "black", linewidth=1, width = 0.4,position= position_dodge(0.9))+
    theme_classic (base_size = 10) + #theme_minimal for all species for mean treatment
    theme (plot.title = element_text ( size = 12), #hjust = 0.5,
           panel.background = element_blank(),
           legend.title = element_text(hjust=0.5),
           legend.key.size = unit(1,"line"),
           #legend.key.spacing = unit(5, "pt"),
           legend.margin = margin(0,0,0,0),
           axis.text=element_text( color = "black"), #angle = 30, hjust = 0.9,
           axis.title= element_blank(), 
           legend.position = "right")
  