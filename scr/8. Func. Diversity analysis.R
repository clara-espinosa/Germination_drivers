library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis):library(purrr);library(stringr):library(tibble)
###################################  MEDITERRANEAN ###############################################################
FD.med %>%
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  dplyr::select(trait, plot, data_type, MPD)-> FD.M.diff

#PAIRED t-test?? http://www.sthda.com/english/wiki/paired-samples-t-test-in-r  

ttest.fd.diff <- function(x) {
  t.test(MPD ~ trait, data = x, paired = TRUE) -> m1 #paired = TRUE, 
  broom::tidy(m1)
}

FD.M.diff%>%
  group_by(data_type)%>%
  do(ttest.fd.diff(.)) %>%
  write.csv("results/FD diff MED.csv")

# with both data types significant differences!
x11()
FD.med%>% # from script 6 coomunity metrics calculations
  filter(!trait == "Germination traits")%>%
  filter (!trait == "Plant traits")%>% 
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cadetblue","cadetblue3","cadetblue1",  
                               "antiquewhite3","antiquewhite1", "cornsilk2","burlywood2", "cornsilk4"))+
  scale_y_continuous(limits = c(0, 1))+
  #coord_flip()+
  facet_grid(~data_type)+
  labs(title = "Functional divergence (MPD)", subtitle = "Mediterranean", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12, angle = 20, vjust = 0.8),
        axis.title.x = element_blank(),
        legend.position = "none")-> FD.med.sep;FD.med.sep

FD.med%>% # from script 6 coomunity metrics calculations
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  filter(data_type=="Abundance")%>%
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  #geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cyan4", "bisque3"))+
  scale_y_continuous(limits = c(0.14, 0.8))+
  geom_segment (aes(x= 1.1,xend =1.9,  y = 0.75, yend= 0.75), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.78, label = "***", size= 4.5)+
  labs()+
  theme_classic(base_size = 16)+
  theme(strip.text = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 10),  #, angle = 20, vjust = 0.8
        axis.text.y = element_text(size =10),
        axis.title= element_blank(),
        legend.position = "none")-> FD.med.ab;FD.med.ab

FD.med%>% # from script 6 coomunity metrics calculations
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  filter(data_type=="Presence")%>%
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  #geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cyan4", "bisque3"))+
  scale_y_continuous(limits = c(0.14, 0.8))+
  geom_segment (aes(x= 1.1,xend =1.9,  y = 0.75, yend= 0.75), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.78, label = "***", size= 4.5)+
  labs()+
  theme_classic(base_size = 16)+
  theme(strip.text = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 10),  #, angle = 20, vjust = 0.8
        axis.text.y = element_text(size =10),
        axis.title= element_blank(),
        legend.position = "none")-> FD.med.pre;FD.med.pre

# combine plots mediterranean
library(patchwork)
FD.med.sep + 
  inset_element(FD.med.ab, left = 0.2, bottom = 0.7, right = 0.494, top = 1)+
  inset_element(FD.med.pre, left = 0.7, bottom = 0.7, right = 1, top = 1)-> graph.FD.med; graph.FD.med


x11()

  

## 5.5.5. explore correlation with plot richness supplementary? ####
#medFD$nbsp from scrp community level analysis
# richness x environmental gradients
as.data.frame(medFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = medFD$nbsp )%>%
  merge(FD.M, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(~micro_variable, scales = "free")+
  labs(title="Mediterranean plots richness across environmental gradients")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_env_M; richness_x_env_M

ggsave(richness_x_env_M, file = "richness x env Med.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 
# richness x FD indices
as.data.frame(medFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = medFD$nbsp )%>%
  merge(FD.M, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~Func.div, scales = "free_x", ncol=4)+
  labs(title="Mediterranean plots richness correlation with FD indices")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_FD_M; richness_x_FD_M

ggsave(richness_x_FD_M, file = "richness x FD Med.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 

# combine plots
library(patchwork)
richness_x_env_M / richness_x_FD_M + 
  plot_layout(heights = c(1,2), guides = "collect") & theme(legend.position = "bottom")-> richness_M;richness_M

ggsave(filename = "Richness_Med.png", plot =richness_M , path = "results/preliminar graphs/richness", 
       device = "png", dpi = 600)


###################################  TEMPERATE  ##################################################################
FD.tem %>%
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  dplyr::select(trait, plot, data_type, MPD)-> FD.T.diff

#PAIRED t-test?? http://www.sthda.com/english/wiki/paired-samples-t-test-in-r  

ttest.fd.diff <- function(x) {
  t.test(MPD ~ trait, data = x, paired = TRUE) -> m1 #paired = TRUE, 
  broom::tidy(m1)
}

FD.T.diff%>%
  group_by(data_type)%>%
  do(ttest.fd.diff(.)) %>%
  write.csv("results/FD diff MED.csv")
# significnat only with presence absence data

FD.tem%>% # from script 6 coomunity metrics calculations
  filter(!trait == "Germination traits")%>%
  filter (!trait == "Plant traits")%>% 
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cadetblue","cadetblue3","cadetblue1",  
                               "antiquewhite3","antiquewhite1", "cornsilk2","burlywood2", "cornsilk4"))+
  scale_y_continuous(limits = c(0, 1))+
  #coord_flip()+
  facet_grid(~data_type)+
  labs(subtitle = "Temperate", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12, angle = 20, vjust = 0.8),
        axis.title.x= element_blank(),
        legend.position = "none")-> FD.tem.sep;FD.tem.sep

FD.tem%>% # from script 6 community metrics calculations
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  filter(data_type=="Abundance")%>%
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  #geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cyan4", "bisque3"))+
  scale_y_continuous(limits = c(0.18, 0.52))+
  geom_segment (aes(x= 1.1,xend =1.9,  y = 0.45, yend= 0.45), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.49, label = "NS", size= 4.5)+
  labs()+
  theme_classic(base_size = 16)+
  theme(strip.text = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 10),  #, angle = 20, vjust = 0.8
        axis.text.y = element_text(size =10),
        axis.title= element_blank(),
        legend.position = "none")-> FD.tem.ab;FD.tem.ab

FD.tem%>% # from script 6 coomunity metrics calculations
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  filter(data_type=="Presence")%>%
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  #geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cyan4", "bisque3"))+
  scale_y_continuous(limits = c(0.18, 0.52))+
  geom_segment (aes(x= 1.1,xend =1.9,  y = 0.45, yend= 0.45), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.47, label = "*", size= 4.5)+
  labs()+
  theme_classic(base_size = 16)+
  theme(strip.text = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 10),  #, angle = 20, vjust = 0.8
        axis.text.y = element_text(size =10),
        axis.title= element_blank(),
        legend.position = "none")-> FD.tem.pre;FD.tem.pre
x11()

# combine plots temperate
library(patchwork)
FD.tem.sep + 
  inset_element(FD.tem.ab, left = 0.2, bottom = 0.7, right = 0.494, top = 1)+
  inset_element(FD.tem.pre, left = 0.7, bottom = 0.7, right = 1, top = 1)-> graph.FD.tem; graph.FD.tem


# Combine both FD plots
graph.FD.med / graph.FD.tem  +
plot_layout(heights = c(1,1), axis = "collect")-> fig4;fig4

ggsave(filename = "fig4.png", plot =fig4, path = "results/Figures", 
       device = "png", dpi = 600)

## 5.5.5. explore correlation with plot richness  ####
# temFD$nbsp from scrp community level analysis
# richnes x environmental gradients
as.data.frame(temFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = temFD$nbsp )%>%
  merge(FD.T, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Rao.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(~micro_variable, scales = "free")+
  labs(title="Temperate plots richness across environmental gradients")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_env_T; richness_x_env_T

ggsave(richness_x_env_T, file = "richness x env Tem.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 
# richness x FD indices
as.data.frame(temFD$nbsp)%>%
  rownames_to_column(var= "plot")%>%
  mutate (richness = temFD$nbsp )%>%
  merge(FD.T, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  gather (Func.div, value, Germ.Rao.ab:Wplant.Mpd.pa)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=richness, x= value, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=richness, x= value),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~Func.div, scales = "free_x", ncol = 4)+
  labs(title="Temperate plots richness correlation with FD indices")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_FD_T; richness_x_FD_T
                                      

ggsave(richness_x_FD_T, file = "richness x FD Tem.png", 
       path = "results/preliminar graphs/richness", scale = 1,dpi = 600) 

library(patchwork)
richness_x_env_T / richness_x_FD_T + 
  plot_layout(heights = c(1,2), guides = "collect") & theme(legend.position = "bottom")-> richness_T;richness_T

ggsave(filename = "Richness_Tem.png", plot =richness_T , path = "results/preliminar graphs/richness", 
       device = "png", dpi = 600)



