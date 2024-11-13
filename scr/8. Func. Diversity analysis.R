library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis):library(purrr);library(stringr):library(tibble)
###################################  MEDITERRANEAN ###############################################################
## germ traits vs plant traits ####
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
##### partial Mediterranean visualization ####
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


## 5.5.5. explore correlation with plot richness supplementary?? ####
# richness x environmental gradients
med_plot %>% #object in script 5 data matrices handling
  merge(read.csv("data/spatial-survey-header-Med.csv"))%>%
  filter(!plot=="A00")%>% # plot with Microlog data NO iButton
  filter(!plot=="B00")%>%# plot with Microlog data NO iButton
  filter(!plot=="C00")%>%# plot with Microlog data NO iButton
  filter(!plot=="D00")%>%# plot with Microlog data NO iButton
  filter(!plot=="B13")%>% # iButtons could no be recovered
  filter(!plot=="B14")%>% # iButtons could no be recovered
  filter(!plot=="B15")%>% # iButtons could no be recovered
  filter(!plot=="C07")%>%
  dplyr::select(site, plot, N_sp)-> plot_richness_M

FD.med%>%
  merge(plot_richness_M, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  dplyr::select(site, plot, N_sp, trait, data_type, MPD, elevation, FDD, GDD, Snow)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=N_sp, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=N_sp, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~micro_variable, scales = "free", ncol = 1)+
  labs(title="Mediterranean plots richness", subtitle = "correlation with microclimatic gradients")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_env_M; richness_x_env_M
 
# richness x FD indices
FD.med%>%
  merge(plot_richness_M, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Rabinalto", "Canada","Solana","Penouta")) %>%
  rename(Snow=Snw)%>%
  dplyr::select(site, plot, N_sp, trait, data_type, MPD, elevation, FDD, GDD, Snow)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=MPD, x= N_sp, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=MPD, x= N_sp),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(trait~data_type, scales = "free")+
  labs(subtitle= "         correlation with FD traits")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text.y.right = element_text(angle = 0),
        axis.title.x = element_text(size =12),
        legend.position = "bottom")-> richness_x_FD_M; richness_x_FD_M


# combine plots
library(patchwork)
richness_x_env_M + richness_x_FD_M + 
  plot_layout(widths = c(1,2), guides = "collect") & theme(legend.position = "bottom")-> richness_M;richness_M

ggsave(filename = "Richness_Med.png", plot =richness_M , path = "results/Supplementary", 
       device = "png", dpi = 600)


###################################  TEMPERATE  ##################################################################
## germ traits vs plant traits ####
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
##### partial temperate visualization ####

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

### final assemblage for visualization ####
# Combine both FD plots
graph.FD.med / graph.FD.tem  +
plot_layout(heights = c(1,1), axis = "collect")-> fig4;fig4

ggsave(filename = "fig4.png", plot =fig4, path = "results/Figures", 
       device = "png", dpi = 600)

## 5.5.5. explore correlation with plot richness  ####
tem_plot%>%
  merge(read.csv("data/spatial-survey-header-Tem.csv", sep = ";"))%>%
  dplyr::select(site, plot, N_sp)-> plot_richness_T
# richnes x environmental gradients
FD.tem%>%
  merge(plot_richness_T, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  dplyr::select(site, plot, N_sp, trait, data_type, MPD, elevation, FDD, GDD, Snow)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=N_sp, x= value_micro, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=N_sp, x= value_micro),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_wrap(~micro_variable, scales = "free", ncol = 1)+
  labs(title="Temperate plots richness", subtitle = "correlation with microclimatic gradients")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")-> richness_x_env_T; richness_x_env_T


# richness x FD indices
FD.tem%>%
  merge(plot_richness_T, by= "plot")%>%
  mutate(site = as.factor(site))%>%
  mutate(site = fct_relevel(site,"Los Cazadores", "Hou Sin Tierri","Los Boches","Hoyo Sin Tierra"))%>% 
  rename(Snow=Snw)%>%
  dplyr::select(site, plot, N_sp, trait, data_type, MPD, elevation, FDD, GDD, Snow)%>%
  gather (micro_variable, value_micro, elevation:Snow)%>%
  ggplot()+
  geom_point(aes(y=MPD, x= N_sp, fill = site), color= "black",shape = 21, size =3)+
  scale_fill_manual(values = c("limegreen","deeppink4","darkorange1", "dodgerblue4"))+
  geom_smooth (aes(y=MPD, x= N_sp),method = "lm", se=T, color = "black", inherit.aes=F)+
  facet_grid(trait~data_type, scales = "free")+
  labs(subtitle= "         correlation with FD traits")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =12),
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text.y.right = element_text(angle = 0),
        axis.title.x = element_text(size =12),
        legend.position = "bottom")-> richness_x_FD_T; richness_x_FD_T
                                    

library(patchwork)

richness_x_env_T + richness_x_FD_T + 
  plot_layout(widths = c(1,2), guides = "collect") & theme(legend.position = "bottom")-> richness_T;richness_T

ggsave(filename = "Richness_Tem.png", plot =richness_T , path = "results/Supplementary", 
       device = "png", dpi = 600)



