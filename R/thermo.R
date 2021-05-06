#load tidyverse library
pacman::p_load(tidyverse, CHNOSZ, cowplot)

#estimate CO activity 
#convert from nM to M and estimate activity from averaged log concentration
CO <- read_csv(paste0(data_path, "geochemistry/geochem.csv")) %>%
  group_by(site) %>%
  summarise(CO = mean(CO, na.rm = T)/1000000000) %>%
  filter(!is.na(site)) %>%
  mutate(log_act = log10(CO),
         species = "CO", 
         molality = CO) %>%
  select(site, species, log_act, molality) #molality here is actually molarity, converted to molality later

#read_activity data - missing CO(aq)
files <- list.files(paste0(data_path, "geochemistry"), full.names = T, pattern = ".*nobalance.txt")

read_files <- function(file){
  site <- str_extract(file, "(?<=data/geochemistry/)(.*)(?=_output_April18_nobalance[.]txt)")
  
  density <-   file %>% 
    read_delim(delim = " ", skip = 10, trim_ws = T, col_names = F) %>%
    filter(X2 == "density") %>%
    select(X4) %>%
    rename(soln_density = "X4") %>%
    mutate(soln_density = as.numeric(soln_density))
  
  
  file %>% 
    read_delim(delim = " ", skip = 41, trim_ws = T, comment = " ----------------------------------------------------------------")  %>%
    mutate(site = site) %>%
    rename(mg_per_kg_soln = "molality") %>%
    rename(molality = "species") %>%
    rename(species = "Aqueous") %>%
    rename(act_coef = `mg/kg`) %>%
    rename(log_act = "sol'n") %>%
    select(site, species, log_act, molality) %>%
    filter(species %in% c("Ca++", "CH3COO-", "CH4(aq)", "Fe++", "H+", "H2(aq)", "HCO3-",
                          "HS-", "Mn++", "NH4+", "NO2-", "NO3-", "SO4--", "O2(aq)",
                          "Cl-", "ClO4-", "ClO3-") & !is.na(log_act)) %>%
    mutate(molality = as.numeric(molality)) %>%
    bind_cols(density)
  
}

file_list = lapply(files, read_files)

activities <- do.call(rbind, file_list) 

densities <- activities %>%
  select(site, soln_density) %>%
  distinct()

activities <- activities %>%
  mutate(log_act = as.numeric(log_act)) %>%
  bind_rows(CO %>% left_join(densities)) %>%
  mutate(molality = if_else(species == "CO", molality/((soln_density*1000 - 28.01*molality)/1000), molality), #convert CO from molarity to molality 
         species = str_remove(species, "[(]aq[)]"),
         species = str_replace(species, "[+][+]", "+2"),
         species = str_replace(species, "--", "-2"),
         species = str_replace(species, "CH3COO-", "acetate")) %>%
  rename(react_prod = "species") %>%
  select(-soln_density)


reactions <- read_csv(paste0(data_path, "geochemistry/reactions.csv"), col_types = cols (.default = "c")) 

#add ferrihydrite to database
ferrihydrite <- mod.OBIGT("ferrihydrite", G=-111200, H=-127800, S=16.7, V=20.88, formula="FeOOH", state="cr", a1.a=8.70, a2.b=36.71, a3.c=-1.0146, E_units = "cal")

thermo_db <- thermo()$OBIGT %>% as_tibble()

#set temperature units to Kelvin
T.units("K")

#set energy units to joules
E.units("J")

logK <- reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)") %>% #pivot from wide to long dataframe
  unite("react_prod", reactant:product, na.rm = TRUE, remove = F) %>% #unite the reactants and products into one column 
  filter(!react_prod == "") %>% #remove any rows with missing react_prod values
  mutate(coeff = as.numeric(coeff)) %>%
  group_by(rxn.number) %>% #group by reaction number for calculations
  summarise(D1 = subcrt(react_prod, coeff, state, T=283.45)$out$logK,#calculate logK using in situ DeMMO1 temperature for all other reactions 
            D2 = subcrt(react_prod, coeff, state, T=285.55)$out$logK,
            D3 = subcrt(react_prod, coeff, state, T=289.35)$out$logK,
            D4 = subcrt(react_prod, coeff, state, T=295.65)$out$logK,
            D5= subcrt(react_prod, coeff, state, T=304.85)$out$logK,
            D6 = subcrt(react_prod, coeff, state, T=294.65)$out$logK) %>%
  pivot_longer(cols = D1:D6, names_to = "site", values_to = "LogK") #pivot from wide to long for merging later
# 
# test2 <- subcrt(species = c("O2", "Fe+2", "H2O", "ferrihydrite", "H+"), 
#                 coeff = c(-1, -4, -6, 4, 8), 
#                 state = c("aq", "aq", "liq", "cr", "aq"), 
#                 T=283.45,
#                 P = 1.013)


logQ_data <- reactions %>%
  pivot_longer(reactant.a:state.i,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)") %>% #pivot from wide to long 
  unite("react_prod", reactant:product, na.rm = TRUE, remove = F) %>% #unite the reactant and product columns into one column called react_prod 
  left_join(activities) %>%
  filter(!is.na(log_act)) %>% #remove any activities with NA values 
  mutate(coeff = as.numeric(coeff),
         logQ = if_else(is.na(reactant), abs(coeff)*log_act, -abs(coeff)*log_act)) %>% #multiply each species by its stoichiometric coefficient 
  group_by(rxn.number,site) %>% #group on the reaction number and site 
  summarise(logQ = sum(logQ), e.donor = unique(e.donor), e.acceptor = unique(e.acceptor)) 


logQ <- logQ_data %>%
  inner_join(activities %>% rename(e.donor = "react_prod", ED_log_act = "log_act", ED_molality = "molality")) %>%
  full_join(logQ_data %>% left_join(activities %>% rename(e.acceptor = "react_prod", EA_log_act = "log_act", EA_molality = "molality")))


deltaG <- logK %>%
  left_join(logQ) %>% #join the logK and logQ tables 
  left_join(reactions %>% select(rxn.number, e.transfer, coeff.a, coeff.b, iron)) %>% #add the reaction number, number of electrons transferred, and minerals from each reaction 
  mutate(e.transfer = as.numeric(e.transfer),
         deltaG = (-2.303*8.314*283.45*(LogK-logQ))/(e.transfer*1000),
         EA = log10(-(deltaG/abs(as.numeric(coeff.a)))*(10^EA_log_act)),
         ED = log10(-(deltaG/abs(as.numeric(coeff.b)))*(10^ED_log_act))) %>%
  filter(!is.na(deltaG)) %>% #calculate deltaG for each reaction at each site
  mutate(EA_availability = abs(as.numeric(coeff.a))*EA_molality,
         ED_availability = abs(as.numeric(coeff.b))*ED_molality,
         energy_density = if_else(EA_availability < ED_availability & !is.na(EA), EA, ED),
         energy_density = if_else(is.na(energy_density) & !is.na(EA), EA, energy_density))  

#write supp table 6
deltaG %>%
  select(rxn.number, site, e.donor, e.acceptor, e.transfer, LogK, logQ, deltaG, energy_density) %>%
  arrange(as.numeric(rxn.number), site) %>%
  write_csv(paste0(write_path, "supp_table6.csv"))

endergonic_rxns <- deltaG %>%
  filter(is.na(energy_density)) %>%
  select(site, rxn.number,e.donor, e.acceptor, iron) %>%
  distinct()
# plot the data! ----------------------------------------------------------

thermo_palette <- c("#f8766d", "#d89000", "#a3a500", "#39b600", "#00bf7d", "#00bfc4", "#00b0f6", "#9590ff", "#e76bf3", "#363534", "brown", "gray", "gray", "gray")
names(thermo_palette) <- c("acetate", "CH4", "ClO4-", "CO", "H2", "HCO3-", "HS-", "NH4+", "NO3-", "O2", "SO4-2", "sulfur", "pyrite", "siderite")

deltaG_plot <- deltaG %>%
  select(rxn.number, site, e.acceptor, e.donor, iron, deltaG, EA, ED)%>%
  pivot_longer(deltaG:ED, names_to = "name", values_to = "value", values_drop_na = T) %>%
  ggplot(aes(value, reorder(rxn.number, -value), shape=site, group=rxn.number)) +
  theme_gray() +
  geom_line(aes(color=iron), size=2.5, alpha=0.6) + #color each line spanning the deltaG values for the six sites by the mineral in the reaction
  geom_point() + 
  scale_shape_manual(values = c(0,1,2,15,16,17)) + #manually set the shapes for each point to denote the six different sites 
  scale_color_manual(values = thermo_palette) +
  #scale_x_reverse() + #reverse the x-axis to show exergonic values on the right, this is standard for this kind of data
  labs(x=expression(Delta~G[r]~'kJ/mol'~e^{textstyle("-")})) + #generate the axis labels
  ylab("Reaction #") +
  geom_vline(xintercept = 0, linetype="dotted", color = "black") + #add a vertical line at zero for reference 
  theme(legend.position = c(.08, .6), legend.text=element_text(size=8), legend.title = element_text(size=8, face="bold")) + #position the legend on the left 
  theme(legend.key.size =  unit(0.1, "in")) +
  facet_wrap(~name, scales = "free_x")#resize the legend to make it fit 

geochemistry <- read_csv(paste0(data_path, "geochemistry/geochem_April2018.csv"))

#plot energy density
energy_density_plot <- deltaG %>%
  select(rxn.number, site, e.acceptor, e.donor, iron, deltaG, energy_density)%>%
  mutate(e_donor_acceptor = if_else(iron %in% c("Fe+2", "pyrite", "siderite"), e.acceptor, e.donor),
         type = if_else(iron %in% c("ferrihydrite", "Fe+2"), "planktonic", "biofilm"),
         iron = if_else(iron != "Fe+2", substr(iron, 1, 3), iron),
         id = paste0(iron, "_", rxn.number)) %>%
  filter(!is.na(energy_density) & e_donor_acceptor != "ClO4-") %>%
  ggplot(aes(energy_density, reorder(id, energy_density), shape=site, color=e_donor_acceptor, label = iron)) +
  theme_gray() +
  geom_point(size = 4, stroke = 1) + 
  scale_shape_manual(values = c(0,1,2,15,16,17)) + #manually set the shapes for each point to denote the six different sites 
  scale_color_manual(values = thermo_palette) +
  labs(x=expression(Delta~G[r]~log~'J/kg'~H[2]~'O')) + #generate the axis labels
  ylab("Reaction #") +
  theme(legend.key.size =  unit(0.1, "in"),
        text=element_text(size=18)) +
  facet_wrap(~type, scales = "free_y")
