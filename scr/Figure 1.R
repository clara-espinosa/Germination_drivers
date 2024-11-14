library (tidyverse); library (dplyr);library(rstatix); library(ggrepel)
library (ggpubr);library (vegan); library (patchwork)
library(cowplot); library(magick); library(png); library(gridExtra)

# figure 1 combining maps with figure with tables
# work in progress!!!
# table with summary of environmental variables
plot_x_env_M2 %>%
  mutate(community = "Mediterranean")%>%
  rbind(plot_x_env_T2 %>%
  mutate(community = "Temperate"))%>%
  group_by(community)%>%
  get_summary_stats()%>%
  dplyr::select(community, variable, mean, max,min,n  )%>%
  gather(statistic, value, mean:n)%>%
  mutate(statistic = as.factor(statistic))%>%
  mutate(statistic = fct_relevel(statistic, "n", "min", "max", "mean"))%>%
  mutate(variable = fct_recode(variable, "Snow" = "Snw", "Elevation"= "elevation"))%>%
  ggplot()+
  geom_text(aes(y= variable, x= statistic, label = round(value, digits=1), vjust=0.5))+
  facet_wrap(~community)+
  #geom_hline(yintercept = 120, linetype = "solid") +
  geom_segment (aes(x= 0.5,xend =4.5,  y = 7.25, yend= 7.25), color = "black", linewidth = 0.8, show.legend = F)+
  annotate ("text", x= c(1), y= 7.5, label = "N plots", size= 4.5)+
  annotate ("text", x= c(2), y= 7.5, label = "Min value", size= 4.5)+
  annotate ("text", x= c(3), y= 7.5, label = "Max value", size= 4.5)+
  annotate ("text", x= c(4), y= 7.5, label = "Mean value", size= 4.5)+
  ggtitle("Environmental variables ranges")+
  theme_classic(base_size= 11)+
  theme(plot.title = element_text(size=14,hjust =-0.06),
        plot.margin = margin (0,0,0,0, unit = "pt"),
    strip.text = element_text(size=14),
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks  = element_blank(),
    axis.title.y  = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.y= element_text(color= "black", size= 12),
    axis.text.x = element_blank())-> table;table

plot_x_env_M2 %>%
  mutate(community = "Mediterranean")%>%
  rbind(plot_x_env_T2 %>%
          mutate(community = "Temperate"))%>%
  group_by(community)%>%
  get_summary_stats()%>%
  dplyr::select(community, variable, mean, max,min,n  )%>%
  gather(statistic, value, mean:n)%>%
  mutate(statistic = as.factor(statistic))%>%
  mutate(statistic = fct_relevel(statistic, "n", "min", "max", "mean"))%>%
  mutate(variable = fct_recode(variable, "Snow" = "Snw", "Elevation"= "elevation"))

# pca figures
bioclim # bioclimatic indices for both communities
germ.trait.pca # all traits pca

# read images (maps from QGis)
EU <- image_read(png::readPNG("images/mapa europa rec.png", native= TRUE))%>%
  image_ggplot()
communities <- image_read(png::readPNG("images/mapa comunidades.png", native= TRUE))%>%
  image_ggplot()
cross <- image_read(png::readPNG("images/mapa cruces.png", native= TRUE))%>%
  image_ggplot()

# combine graphs
EU+communities+cross-> maps;maps
ggsave(filename = "Figure 1A.png", plot =maps , path = "results/Figures", 
       device = "png", dpi = 600) 
  
maps/bioclim / germ.trait.pca  +  plot_layout(heights =  c(1,1,1))

cowplot::plot_grid(bioclim ,germ.trait.pca,  ncol = 1)

