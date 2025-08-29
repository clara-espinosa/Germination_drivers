library (binom);library (ggsignif);library (rstatix); library (stringr);
library(patchwork);library(scales);library (tidyverse)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")
detach(package:plyr)

### MACROCLIMATE ####
read.csv("data/wp_villa_2020_2025.csv", sep = ",") %>%
  #mutate(Time = strptime(as.character(Time), "%Y/%m/%d %H:%M"))%>% #specify format of Time variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>% 
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  dplyr::filter(Micro == "central") %>% # keep logger from center of the plot
  na.omit()%>%
  filter(Year > 2020) %>% 
  mutate (Site_year = paste (Site, Year)) %>%
  mutate(WP = rowMeans((cbind(wp1, wp2)))) %>% # calculate mean value between 2-3 gypsum sensors
  dplyr::select(! c(wp1, wp2)) %>% 
  group_by(Community, Site, Month, Day) %>% # calculate mean values x day
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time), 
            WPmean = mean(WP), WPmax = max(WP), WPmin = max (WP), abs_WP = sum(WP)) %>%
  mutate(GDD = ifelse(T >= 5, T, 0)) %>%
  group_by(Community, Site,  Month) %>% # calculate mean values x month 
  summarise(T = mean(T), X = max(X), N = min(N), abs_WP = sum(abs_WP), GDD = sum(GDD),n = length(n),
            WPmean = mean(WPmean), WPmax = mean(WPmax), WPmin = mean(WPmin))%>%
  print(n=48)-> villa_clima

x11()
villa_clima %>%
  group_by(Month)%>%
  summarise_at(vars(T:WPmin), mean, na.rm = TRUE)%>%
  ggplot() +
  geom_col(aes (x= Month, y=WPmax*2), colour = "black", fill = "darkgoldenrod1")+
  geom_line (aes (x=Month, y=N), colour = "darkgoldenrod", linewidth =0.75) + 
  geom_line (aes (x=Month, y=X), colour = "darkgoldenrod", linewidth =0.75) + 
  geom_ribbon (aes (x=Month, ymin =N, ymax=X), fill = "darkgoldenrod", alpha =0.3) + 
  scale_x_continuous (limits = c(0.5,12.5), breaks = seq (1, 12, by= 1))+
  scale_y_continuous(limits = c(-10, 45), breaks = seq (-10, 40, by = 10),
                     sec.axis = sec_axis(transform = ~./20, name = expression(paste("Mean Maximum Water Potential (MPa)")))) + 
  geom_hline(yintercept=30, linetype ="dashed", linewidth =0.5, colour = "darkred")+
  ggthemes::theme_tufte(base_size = 6) + 
  theme (plot.margin = margin(0, 0, 0, 0, "pt"),
         text = element_text(family = "sans"),
         panel.background = element_rect(color = "black", fill = NULL),
         plot.title = element_text (size = 6, margin=margin(0,0,0,0), color = "darkgoldenrod1", face= "bold"), #hjust = 0.5,
         axis.title.y = element_text (size=5), 
         axis.text.y.right = element_blank(),
         axis.ticks.y.right=element_blank(),
         legend.title = element_text(size =5),
         legend.text = element_text (size =4)) +
  labs (title = "Mediterranean Alpine", y= "Temperature (ºC)", x = "Month")->fig1a;fig1a

### TEMPERATE CLIMATE ####
read.csv("data/wp_picos_2021_2025.csv", sep = ",") %>%
  #mutate(Time = strptime(as.character(Time), "%Y/%m/%d %H:%M"))%>% #specify format of Time variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>% 
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  dplyr::filter(Micro == "central") %>% # keep logger from center of the plot
  na.omit()%>%
  filter(!Site == "El Cable")%>%
  filter(!Site == "Urriellu")%>%
  filter(!Site == "Cabaña Verónica")%>%
  filter(!Site == "Horcados Rojos")%>%
  filter(Year > 2020) %>% 
  mutate (Site_year = paste (Site, Year)) %>%
  mutate(WP = rowMeans((cbind(wp1, wp2)))) %>% # calculate mean value between 2-3 gypsum sensors
  dplyr::select(! c(wp1, wp2, wp3)) %>% 
  group_by(Community, Site, Month, Day) %>% # calculate mean values x day
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time), 
            WPmean = mean(WP), WPmax = max(WP), WPmin = max (WP), abs_WP = sum(WP)) %>%
  mutate(GDD = ifelse(T >= 5, T, 0)) %>%
  group_by(Community, Site,  Month) %>% # calculate mean values x month 
  summarise(T = mean(T), X = max(X), N = min(N), abs_WP = sum(abs_WP), GDD = sum(GDD),n = length(n),
            WPmean = mean(WPmean), WPmax = mean(WPmax), WPmin = mean(WPmin))%>%
  print(n=48)-> picos_clima

x11()
picos_clima %>%
  group_by(Month)%>%
  summarise_at(vars(T:WPmin), mean, na.rm = TRUE)%>%
  ggplot() +
  geom_col(aes (x= Month, y=WPmax*2), colour = "black", fill = "forestgreen")+
  geom_line (aes (x=Month, y=N), colour = "darkgreen", linewidth =0.75) + 
  geom_line (aes (x=Month, y=X), colour = "darkgreen", linewidth =0.75) + 
  geom_ribbon (aes (x=Month, ymin =N, ymax=X), fill = "darkgreen", alpha =0.3) + 
  scale_x_continuous (limits = c(0.5,12.5), breaks = seq (1, 12, by= 1))+
  scale_y_continuous(limits = c(-10, 45), breaks = seq (-10, 40, by = 10),
                     sec.axis = sec_axis(transform = ~./20, name = expression(paste("Mean Maximum Water Potential (MPa)")))) + 
  geom_hline(yintercept=30, linetype ="dashed", linewidth =0.5, colour = "darkred")+
  ggthemes::theme_tufte(base_size = 6) + 
  theme (plot.margin = margin(0, 0, 0, 0, "pt"),
         text = element_text(family = "sans"),
         panel.background = element_rect(color = "black", fill = NULL),
         plot.title = element_text (size = 6, margin=margin(0,0,0,0), color = "forestgreen", face= "bold"), #hjust = 0.5,
         axis.title.y = element_text (size=5), 
         axis.text.y.left = element_blank(),
         axis.ticks.y.left=element_blank(),
         legend.title = element_text(size =5),
         legend.text = element_text (size =4)) +
  labs (title = "Temperate Alpine", y= "Temperature (ºC)", x = "Month")->fig1b;fig1b

# combine panels ###
library(patchwork)
fig1a+ fig1b + 
  plot_layout(widths = c(1,1), axis_titles = "collect")-> fig1;fig1

ggsave(filename = "Fig1clima.png", plot =fig1 , path = "results/figures", 
       device = "png", dpi = 600, width = 70, height = 50, units = "mm") 
