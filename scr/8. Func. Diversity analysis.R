library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis):library(purrr);library(stringr):library(tibble)
###################################  MEDITERRANEAN ###############################################################
FD.med%>% # from script 6 coomunity metrics calculations
  filter(!trait == "Germination traits")%>%
  filter (!trait == "Plant traits")%>% 
  ggplot()+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  #scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  scale_y_continuous(limits = c(0, 1))+
  #coord_flip()+
  facet_grid(~data_type)+
  labs(title = "Mediterranean community FD differences", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12, angle = 20, vjust = 0.8),
        legend.position = "none")-> Full.FD.M;Full.FD.M

# combine plots
library(patchwork)
Full.FD.M/ Full.FD.T+ 
  plot_layout(heights = c(1,1), guides = "collect")->Full.FD; Full.FD

ggsave(filename = "Full FD subsets.png", plot =Full.FD, path = "results/preliminar graphs", 
       device = "png", dpi = 600)
x11()

## 5.5.3 MPD differences between germination and adult plant traits ####
strip_names <- c("ab"= "Abundance data", "pa"= "Presence/absence data", 
                 "Mpd" = "Mean Pairwise Dissimilarity")
FD.M%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", names = c("Traits","indices", "data_type"), too_many = "merge")%>%
  group_by(Traits, indices, data_type)%>%
  get_summary_stats(value)%>%
  mutate(Traits= as.factor(Traits))%>%
  mutate(Traits = fct_recode(Traits, "Germination"="Germ" , "Adult Plant"="Wplant"))%>%
  ggplot()+
  geom_errorbar(aes(Traits, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_point(aes(x= Traits, y= mean, fill = Traits), color = "black", shape= 21, size = 5)+
  scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  facet_grid(~data_type, labeller = as_labeller(strip_names))+
  labs(title = "Mediterranean community FD differences", y= "Mean dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")->FD.traits.M;FD.traits.M

  ggsave(FD.traits.M, file = "FD traits diff Med.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 
  
FD.M %>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", 
                       names = c("Traits","indices", "data_type"), too_many = "merge")-> FD.M.diff

#PAIRED t-test?? http://www.sthda.com/english/wiki/paired-samples-t-test-in-r  

ttest.fd.diff <- function(x) {
  t.test(value ~ Traits, data = x) -> m1 #paired = TRUE, 
  broom::tidy(m1)
}

FD.M.diff%>%
  group_by(indices, data_type)%>%
  do(ttest.fd.diff(.)) %>%
  write.csv("results/FD diff MED.csv")

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
FD.tem%>% # from script 6 coomunity metrics calculations
  filter(!trait == "Germination traits")%>%
  filter (!trait == "Plant traits")%>% 
  ggplot()+
  geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_point(aes(x= trait, y= mean, fill = trait), color = "black", shape= 21, size = 5)+
  #scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  scale_y_continuous(limits = c(0.05, 0.7))+
  coord_flip()+
  facet_grid(~data_type)+
  labs(title = "Temperate community FD differences", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        #axis.title.x = element_blank(),
        legend.position = "none")-> Full.FD.T;Full.FD.T
x11()
## 5.5.3  MPD differences between germination and adult plant traits ####
strip_names <- c("ab"= "Abundance data", "pa"= "Presence/absence data", 
                 "Mpd" = "Mean Pairwise Dissimilarity")
FD.T%>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", names = c("Traits","indices", "data_type"), too_many = "merge")%>%
  group_by(Traits, indices, data_type)%>%
  get_summary_stats(value)%>%
  mutate(Traits= as.factor(Traits))%>%
  mutate(Traits = fct_recode(Traits, "Germination"="Germ" , "Adult Plant"="Wplant"))%>%
  ggplot()+
  geom_errorbar(aes(Traits, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth = 0.5, size =1) +
  geom_point(aes(x= Traits, y= mean, fill = Traits), color = "black", shape= 21, size = 5)+
  scale_fill_manual(values = c("dodgerblue4","limegreen"))+
  facet_grid(indices~data_type, labeller = as_labeller(strip_names), scales = "free_y")+
  labs(title = "Temperate community FD differences", y= "Mean dissimilarity")+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size =16),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_blank(),
        legend.position = "bottom")->FD.traits.T;FD.traits.T

ggsave(FD.traits.T, file = "FD traits diff Tem.png", 
       path = "results/preliminar graphs", scale = 1,dpi = 600) 

FD.T %>%
  gather (Func.div, value, Germ.Mpd.ab:Wplant.Mpd.pa)%>%
  separate_wider_delim(Func.div, delim = ".", 
                       names = c("Traits","indices", "data_type"), too_many = "merge")-> FD.T.diff

ttest.fd.diff <- function(x) {
  t.test(value ~ Traits, data = x) -> m1
  broom::tidy(m1)
}

FD.T.diff%>%
  group_by(indices, data_type)%>%
  do(ttest.fd.diff(.)) %>%
  write.csv("results/FD diff TEM.csv")

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



