library(tidyverse);library(zoo);library (rstatix);library(lubridate)
Sys.setlocale("LC_ALL", "English")

### GERM DRIVERS Mean snow and water potential per year #####
read.csv("data/wp_picos_2016_2021.csv")%>% # data from picos (2016-2021)
  rbind(read.csv("data/wp_picos_2021_2024.csv"))%>% # data from picos (2021-2024)
  rbind(read.csv("data/wp_villa_2020_2024.csv"))%>% # data from villabandin/Omaña (2020-2024)
  select(!wp3)%>%
  #filter(Micro =="central")%>%
  filter(!Site == "El Cable")%>%
  filter(!Site == "Urriellu")%>%
  filter(!Site == "Cabaña Verónica")%>%
  filter(!Site == "Horcados Rojos")%>%
  mutate(Community = as.factor(Community))%>%
  mutate (Community = recode (Community, "Picos"="Temperate" , "Omana"= "Mediterranean" ))%>%
  mutate(Community = fct_relevel (Community, "Temperate","Mediterranean")) %>%
  #mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  mutate(Year = lubridate::year(Time))%>%
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(WP = rowMeans((cbind(wp1, wp2)))) %>% 
  mutate(Temperature = as.numeric(Temperature))%>%
  #filter(Year>2020)%>%
  group_by(Community,Site, Micro, Year, Month, Day) %>%
  summarise(Tmean = mean(Temperature, na.rm = T), Tmax = max(Temperature, na.rm = T), Tmin = min(Temperature, na.rm = T), # Daily mean, max, min
            WPmean = mean(WP, na.rm = T), WPmax = max(WP), WPmin = min (WP), cumWP= sum(WP))%>%
  mutate(Snow = ifelse(Tmax < 0.5 & Tmin > -0.5, 1, 0)) %>% # Day with or without snow
  group_by(Community, Site, Micro, Year, Month) %>%
  summarise(Tmean = mean(Tmean), Tmax = mean(Tmax), Tmin = mean(Tmin), # Daily mean, max, min
            WPmean = mean(WPmean), WPmax = mean(WPmax), WPmin = mean(WPmin), cumWP= sum(cumWP), 
            Snow = sum(Snow), n_days= n_distinct(Day))%>%
  filter(n_days>27)%>% # remove months with too few days with data
  group_by(Community, Site,Micro, Year)%>%
  summarise(Tmean = mean(Tmean), Tmax = mean(Tmax), Tmin = mean(Tmin), # Daily mean, max, min
            WPmean = mean(WPmean), WPmax = mean(WPmax), WPmin = mean(WPmin), cumWP= sum(cumWP), 
            Snow = sum(Snow), n_months=n_distinct(Month))%>%
  filter(n_months >10)%>%
  #print(n=30)
  group_by (Community, Micro, Site)%>%
  summarise (Snow = mean(Snow), cumWP = mean(cumWP))

# GROWING SEASON AND WATER STRESS #####
read.csv("data/wp_picos_2016_2021.csv")%>% # data from picos (2016-2021)
  rbind(read.csv("data/wp_picos_2021_2024.csv"))%>% # data from picos (2021-2024)
  rbind(read.csv("data/wp_villa_2020_2024.csv"))%>% # data from villabandin/Omaña (2020-2024)
  select(!wp3)%>% # remove wp3 (third sensor of water potential from 2021 was not put in the field, results in NAs)
  mutate(Year = lubridate::year(Time)) %>% # create separate columns for Year, Month, Day and Hour to facilitate filtering below
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  mutate (code = paste (Site, Micro, Year))%>% #unique id variable to group our data per site, microclimatic conditions and year
  convert_as_factor(Community, File_name, Site, Micro, ID, code)-> df
str(df)
### growing season determination ####
df %>%  
  mutate(WP = rowMeans((cbind(wp1, wp2)))) %>% # mean value between WP sensors
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  group_by(Community, Site, Micro, Year, code, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature)) %>% # Daily mean
  mutate (data_start = first(Day), data_end = last(Day)) %>% # get starting and ending data points x year
  mutate(t5 = ifelse(T>=5, 1, 0))%>% # days with Tmean>= 5 (Körner 2021) =1
  mutate(length = rollsumr(t5, k = 3, fill= 0)) %>% #sum the 2 previous rows of t5 k=3 porque la literature solo dice "varios"
  filter(! (length < 3)) %>% # Filter date  with 2 consecutive days with Tmean>5ºC (Körner limit) quita las fechas que no son growing season
  group_by(Community, Site, Micro, Year,code) %>% #separate x year
  summarise(GS_start = first (Day), GS_end = last(Day), # get the first and last day of the growing season (puede haber dias entremedias que la temperatura baja de 5 grados investigar)
            data_start = first(data_start), data_end = last(data_end )) %>% # para mirar cuando empieza y termina de coger los datos cada año en cada sitio
  mutate (GS_length = difftime(GS_end , GS_start, units = "days")) %>% # count how many days of growing season
  mutate (data_days = data_end - data_start) %>% # count how many days of data
  select(Community,Site, Micro, Year, code, data_start, data_end, data_days, GS_start, GS_end, GS_length)-> grow

df %>%
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  merge(grow) %>% # first date with more than 2 days with Tmean>5ºC
  filter (! (data_days <300)) %>% # 300 to include data from PICOS in 2019 (check calendar in doc folder)
  filter(Time > (GS_start)) %>% # these two next lines filter for growing season, if you don't want to filter just deactivate
  filter(Time < (GS_end)) %>%
  mutate(WP = rowMeans((cbind(wp1, wp2)))) %>% 
  select(! c(wp1, wp2)) %>% 
  group_by(Community, Site, Micro, Year, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), #n = length(Time), 
            WPmean = mean(WP), WPmax = max(WP), WPmin = min (WP), WPtotal = sum(WP)) %>%
  group_by(Community, Site, Micro) %>% #separate x year 
  summarise(T = mean(T), WPmean = mean(WPmean), WPmax = max(WPmax ), WPmin = min (WPmin), WPtotal = sum(WPtotal))
