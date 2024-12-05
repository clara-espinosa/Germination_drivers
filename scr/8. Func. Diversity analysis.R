library(tidyverse);library(ggplot2);library(dplyr);library(rstatix)
library(lubridate);library(FD); library(psych);library(picante);library(vegan)
library (gawdis):library(purrr);library(stringr):library(tibble)
###################################  MEDITERRANEAN ################################################
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
read.csv("results/GLMs auto FD Med.csv", sep=";") %>% # table with glm model results from script 7
  convert_as_factor(trait, data_type,term) %>%
  filter(!term == "(Intercept)")%>%
  filter(!term == "scale(auto)")%>%
  filter(trait == "germtraits"|trait == "Wplanttraits")%>%
  mutate(trait= fct_recode(trait, "Germination traits"="germtraits","Plant traits"="Wplanttraits"))%>%
  mutate (trait = fct_relevel(trait, "Germination traits","Plant traits"))%>%
  mutate(term= fct_recode(term, "Elevation" = "scale(elevation)", "FDD"="scale(FDD)",
                          "GDD"="scale(GDD)", "Snow"="scale(Snw)"))%>%
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  mutate(CI = 1.96*std.error)%>% # confidence interval multiply 1.96 per std error
  mutate(CImin = estimate-CI, 
         CImax= estimate+CI)%>%
  mutate(alpha = ifelse(p.value<0.05, 1, 0.9))%>%
  ggplot(aes(x= trait, y =estimate, ymin = CImin, ymax = CImax))+
  labs(title = "Functional divergence (MPD)\nMediterranean")+
  geom_errorbar (aes(color=trait, alpha = alpha),width = 0.15, linewidth =1.2, show.legend = F) + #, color="black" 
  geom_point(aes(color=trait, alpha = alpha),size = 3, show.legend = T) +#
  guides(alpha="none")+
  scale_color_manual(values = c("cyan4", "bisque3") )+
  facet_grid (term~data_type,  scales = "free") + #ncol = 8, nrow =2,
  geom_hline(yintercept = 0, linetype = "dashed", linewidth =1, color = "red") +
  coord_flip() +
  labs(y = "GLMs parameter estimate") + # Temperate   Mediterranean
  theme_classic (base_size = 12)+
  #ggthemes::theme_tufte(base_size = 14) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size = 16),
        plot.margin = margin(0,0,0.15,0,unit = "cm"),
        legend.title = element_text(size = 12),
        legend.position = "none",
        strip.background = element_rect(color = "black", fill = NULL),
        strip.text.x = element_text(size =15),
        strip.text.y = element_text(angle = 0, size = 15),
        panel.background = element_blank (),
        plot.tag.position = c(0.015,1),
        #axis.title.x = element_text(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank (),
        axis.text.x = element_text(size = 10, color = "black"),
        #axis.text.y = element_text(size = 10, color = "black"),
        axis.text.y = element_blank (),
        axis.ticks = element_blank())->fd.med.tog.micro;fd.med.tog.micro

FD.med%>% # from script 6 coomunity metrics calculations
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  #filter(data_type=="Abundance")%>%
  ggplot()+
  labs(y= "Mean pairwise dissimilarity")+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), show.legend = F, color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(name = "Traits",labels = c("Germination\ntraits", "Plant traits"),values = c("cyan4", "bisque3"))+
  scale_y_continuous(limits = c(0.14, 0.8))+
  facet_grid(~data_type)+
  geom_segment (aes(x= 1.1,xend =1.9,  y = 0.75, yend= 0.75), color = "black", linewidth = 1, show.legend = F)+
  annotate ("text", x= 1.5, y= 0.78, label = "***", size= 4.5)+
  labs()+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size= 15),
        plot.subtitle = element_text (hjust=0.5),
        plot.margin = margin(0,0.15,0.5,0, unit = "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12),  #, angle = 20, vjust = 0.8
        axis.text.y = element_text(size =12),
        axis.title.y= element_text(size =14),
        axis.title.x = element_blank(),
        legend.position.inside = c(1.1, 0.5), 
        legend.position = "none")-> FD.med.tog;FD.med.tog

#combine plots
library(patchwork)
fd.med.tog.micro/FD.med.tog + plot_layout(heights = c(0.25,0.75), guides= "collect")


## supplementary figure S2?#####
FD.med%>% # from script 6 coomunity metrics calculations
  filter(!trait == "Germination traits")%>%
  filter (!trait == "Plant traits")%>% 
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  ggplot()+
  labs(title = "Functional divergence (MPD)\nMediterranean", y= "Mean pairwise dissimilarity")+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 1, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cadetblue","cadetblue3","cadetblue1",  
                               "antiquewhite3","antiquewhite1", "cornsilk2","burlywood2", "cornsilk4"))+
  scale_y_continuous(limits = c(0, 1))+
  #coord_flip()+
  facet_grid(~data_type)+
  #labs(title = "Functional divergence (MPD)", subtitle = "Mediterranean", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 14)+
  theme(strip.text = element_text(size=15),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12, angle = 20, vjust = 0.8),
        axis.title.x = element_blank(),
        legend.position = "none")-> FD.med.sep;FD.med.sep



# combine plots mediterranean
library(patchwork)

# version 2
(FD.med.tog /FD.med.sep )+ plot_layout(heights = c(1,0.5))-> graph.FD.med; graph.FD.med

#version 3

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


###################################  TEMPERATE  ##################################################
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
read.csv("results/GLMs auto FD Tem.csv", sep=",") %>% # table with glm model results from script 7
  convert_as_factor(trait, data_type,term) %>%
  filter(!term == "(Intercept)")%>%
  filter(!term == "scale(auto)")%>%
  filter(trait == "germtraits"|trait == "Wplanttraits")%>%
  mutate(trait= fct_recode(trait, "Germination traits"="germtraits","Plant traits"="Wplanttraits"))%>%
  mutate (trait = fct_relevel(trait, "Germination traits","Plant traits"))%>%
  mutate(term= fct_recode(term, "Elevation" = "scale(elevation)", "FDD"="scale(FDD)",
                          "GDD"="scale(GDD)", "Snow"="scale(Snw)"))%>%
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  mutate(CI = 1.96*std.error)%>% # confidence interval multiply 1.96 per std error
  mutate(CImin = estimate-CI, 
         CImax= estimate+CI)%>%
  mutate(alpha = ifelse(p.value<0.05, 1, 0.9))%>%
  ggplot(aes(x= trait, y =estimate, ymin = CImin, ymax = CImax))+
  geom_errorbar (aes(color=trait, alpha = alpha),width = 0.15, linewidth =1.2) + #, color="black" 
  geom_point(aes(color=trait, alpha = alpha),size = 3, show.legend = F) +#
  guides(alpha="none")+
  scale_color_manual(values = c("cyan4", "bisque3") )+
  facet_grid (term~data_type,  scales = "free") + #ncol = 8, nrow =2,
  geom_hline(yintercept = 0, linetype = "dashed", linewidth =1, color = "red") +
  coord_flip() +
  labs(title = "Temperate",y = "GLMs parameter estimate" )+
  theme_classic (base_size = 12)+
  #ggthemes::theme_tufte(base_size = 14) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size = 16),
        plot.margin = margin(0,0,0.15,0,unit = "cm"),
        legend.title = element_text(size = 12),
        legend.position = "none",
        strip.background = element_rect(color = "black", fill = NULL),
        strip.text.x = element_text(size =15),
        strip.text.y = element_text(angle = 0, size = 15),
        panel.background = element_blank (),
        plot.tag.position = c(0.015,1),
        #axis.title.x = element_text(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank (),
        axis.text.x = element_text(size = 10, color = "black"),
        #axis.text.y = element_text(size = 10, color = "black"),
        axis.text.y = element_blank (),
        axis.ticks = element_blank())->fd.tem.tog.micro;fd.tem.tog.micro

# significances for grouped traits

sig.segment <- as.data.frame(cbind(data_type = c("Presence/Absence data", "Abundance data"),
                                   x= c(1.1,1.1), xend=c(1.9,1.9), y= c(0.45, 0.45), yend=c(0.45, 0.45)))
sig.segment%>%
  mutate(x= as.numeric(x),
         xend= as.numeric(xend),
         y= as.numeric(y),
         yend= as.numeric(yend))->sig.segment
str(sig.segment)

sig.label <- as.data.frame(cbind(data_type = c("Presence/Absence data", "Abundance data"),
                                 x= c(1.5, 1.5), y= c(0.47, 0.47), label = c("*", "NS")))
sig.label%>%
  mutate(x= as.numeric(x),
         y= as.numeric(y))->sig.label
str(sig.label)

FD.tem%>% # from script 6 coomunity metrics calculations
  filter(trait == "Germination traits"|trait == "Plant traits" )%>%
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  #filter(data_type=="Abundance")%>%
  ggplot()+
  labs(y= "Mean pairwise dissimilarity")+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), show.legend = F, color = "black", shape= 21, size = 5, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cyan4", "bisque3"))+
  scale_y_continuous(limits = c(0.18, 0.52))+
  facet_grid(~data_type)+
  geom_segment (data=sig.segment, aes(x= x,xend =xend,  y = y, yend= yend), color = "black", linewidth = 1, show.legend = F)+
  geom_text (data= sig.label, aes(x= x, y= y, label = label), size= 4.5)+
  labs()+
  theme_classic(base_size = 16)+
  theme(strip.text = element_text(size= 15),
        plot.subtitle = element_text (hjust=0.5),
        plot.margin = margin(0,0.15,0,0, unit = "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12),  #, angle = 20, vjust = 0.8
        axis.text.y = element_text(size =12),
        axis.title.y= element_text(size =14),
        axis.title.x = element_blank(),
        legend.position = "none")-> FD.tem.tog;FD.tem.tog
#combine plots
library(patchwork)
fd.tem.tog.micro/FD.tem.tog + plot_layout(heights = c(0.25,0.75))-> fd.tem.pooled;fd.tem.pooled

# supplementary figure s2 ####

FD.tem%>% # from script 6 coomunity metrics calculations
  filter(!trait == "Germination traits")%>%
  filter (!trait == "Plant traits")%>% 
  mutate(data_type = fct_recode(data_type,"Presence/Absence data"= "Presence","Abundance data"= "Abundance"))%>%
  ggplot()+
  labs(title = "Temperate", y= "Mean pairwise dissimilarity")+
  #geom_errorbar(aes(trait, mean, ymin = mean-ci, ymax = mean+ci), color = "black",linewidth =1) +
  geom_boxplot(aes(x= trait, y= MPD, fill = trait), color = "black")+
  geom_point(aes(x= trait, y= MPD, fill = trait), color = "black", shape= 21, size = 1, alpha = 0.5,position = "jitter")+
  scale_fill_manual(values = c("cadetblue","cadetblue3","cadetblue1",  
                               "antiquewhite3","antiquewhite1", "cornsilk2","burlywood2", "cornsilk4"))+
  scale_y_continuous(limits = c(0, 1))+
  #coord_flip()+
  facet_grid(~data_type)+
  #labs(subtitle = "Temperate", y= "Mean pairwise dissimilarity")+
  theme_classic(base_size = 14)+
  theme(strip.text = element_blank(),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size= 12, angle = 20, vjust = 0.8),
        axis.title= element_text(),
        legend.position = "none")-> FD.tem.sep;FD.tem.sep


x11()

# combine plots temperate
library(patchwork)

# version 2
(FD.tem.tog/FD.tem.sep)+ plot_layout(heights = c(1,0.5))-> graph.FD.tem; graph.FD.tem

# combine supplementary figure s2
(FD.med.sep/FD.tem.sep)+ plot_layout(heights = c(1,1))-> graph.FD.each;graph.FD.each
ggsave(filename = "Figure S2.png", plot =graph.FD.each, path = "results/Supplementary", 
       device = "png", dpi = 600)

### final assemblage for visualization ####
# Combine both FD plots
graph.FD.med / graph.FD.tem  +
plot_layout(heights = c(.75,0.25,1), axis = "collect")

# version 3
fd.med.tog.micro/FD.med.tog/ fd.tem.tog.micro/FD.tem.tog +
  plot_layout (heights = c(0.2,0.3,0.2,0.3))-> fig4;fig4
  

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




