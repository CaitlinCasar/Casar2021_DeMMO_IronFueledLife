pacman::p_load(tidyverse, lubridate, zoo)

#measurement units:
#flow rate mL/min
#temp C
#cond uS
#TDS ppm
#gas nM
#DIC mM 
#ions/metals mg/L except S2 which is ug/L

geochemistry <- read_csv("../data/geochemistry/geochem.csv") %>%
  mutate(date = as.Date(mdy(date), "%m-%d-%Y"),
         date =as.yearmon(date, "%m/%Y"),
         S2 = S2*0.001) %>% #convert units to mg/L
  filter(date == "April 2018") %>%
  #select(site, Ca, Mg, Na, Fe, Mn, Si, Cl)
  select(site, temp, ORP, DOC, NO3, NH4, Fe2, S2, DO, SO4, Cl)

#average gas measurements and convert from nM to mg/L
averaged_data <- read_csv("../data/geochemistry/geochem.csv") %>%
  mutate(CH4 = (CH4*10^-9)*16.04*1000,
         H2 = (H2*10^-9)*1.00794*1000,
         CO2 = (CO2*10^-9)*44.01*1000,
         CO = (CO*10^-9)*28.01*1000) %>%
  group_by(site) %>%
  summarise(CH4 = mean(CH4, na.rm = T),
            H2 = mean(H2, na.rm = T),
            CO2 = mean(CO2, na.rm = T),
            CO = mean(CO, na.rm = T),
            pH = mean(pH, na.rm = T))
  

geochemistry <- geochemistry %>%
  left_join(averaged_data)

write_csv(geochemistry, "../data/geochemistry/geochem_April2018.csv")

#plot geochem 
geochem_plot <- geochemistry %>%
  pivot_longer(-site, names_to = "name", values_to = "value") %>%
  mutate(name = factor(name, levels=c("H2", "CO", "CH4", "S2", "NH4", "NO3", "DOC", "Fe2", "SO4"))) %>%
  filter(!is.na(name)) %>%
  ggplot(aes(name, value, color=name)) +
  geom_point(size=5) + 
  scale_y_log10() + 
  #ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-std_dev, ymax=mean+std_dev), width=0.2) + 
  ggplot2::coord_flip() +
  ggplot2::geom_line(ggplot2::aes(group=name), color="gray", linetype="dotted") +
  facet_wrap(~site) +
  theme(legend.position = "none", 
        axis.title.y=ggplot2::element_blank()) +
  ylab("Concentration (mg/L)") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) 
