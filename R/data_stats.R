pacman::p_load(tidyverse, cowplot)

# calculate metagenome stats ----------------------------------------------------

metagenome_data <- read_csv("../data/metagenome_data.csv") %>%
  mutate(check = if_else(call == "ambiguous" & str_detect(category, "Metal reduction") & str_detect(category, "iron_reduction"), "fegenie", call),
         category = if_else(call == "ambiguous" & str_detect(category, "Metal reduction"), str_split(category, ",")[[1]][1], category),
         call = check) %>% #re-annotate mismatching annotations as "ambiguous
  select(-check) %>%
  group_by(site, genome, call, category, gene_function, Gene.abbreviation, Gene.name) %>%
  summarise(rel_hits = sum(rel_hits), hits = sum(hits))

gene_counts <- read_delim("../data/metabolic/metagenomes/geneCounts.txt", delim = "\t", col_names = F) %>%
  separate(X1, c("id", "gene_counts"), sep = " ") %>%
  separate(id, c("site", "genome")) %>%
  mutate(site = str_remove(site, "eMMO"))

annotations <- metagenome_data %>%
  group_by(call, site) %>%
  summarise(hits = sum(hits)) %>%
  group_by(site) %>%
  mutate(total_hits = sum(hits)) %>%
  left_join(gene_counts) %>%
  group_by(call, site) %>%
  mutate(prop_hits = hits/total_hits,
         rel_hits = hits/as.numeric(gene_counts))

annotations_plot <- annotations %>%
  pivot_longer(prop_hits:rel_hits, names_to = "name", values_to = "value") %>%
  mutate(name = str_to_title(str_replace_all(name, "_", " ")),
         site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")))) %>%
  ggplot(aes(site, value, fill = call)) +
  geom_bar(stat="identity") +
  facet_wrap(~name, scales = "free") + 
  coord_flip() +
  ggtitle("A. Metagenome annotations") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        text = element_text(size = 18))

metagenome_stats_table <- gene_counts %>%
  mutate(site = factor(site, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC"))) %>%
  select(-genome) %>%
  pivot_longer(-site, names_to = "observation", values_to = "value") %>%
  mutate(observation = str_to_title(str_replace_all(observation, "_", " "))) %>%
  ggplot(aes(x=observation,y=site, label = value)) + 
  geom_tile(fill = "white") + geom_text(colour = "black") +
  facet_grid(vars(site),vars(observation), scales = "free", labeller = label_wrap_gen(width = 2, multi_line = TRUE)) +
  ggtitle("C. Metagenome Stats") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        axis.ticks=element_blank(),
        text = element_text(size = 18))

# calculate MAG stats ----------------------------------------------------
files <- list.files("../data/metabolic/genomes", full.names = T, pattern = ".*checkm")
read_files <- function(file){
  file %>%
    read_delim(delim = "\t")
}
checkm_list = lapply(files, read_files)
checkm_stats <- reduce(checkm_list, full_join) %>%
  separate(`Bin Id`, into = c("site", "genome")) %>%
  mutate(site = str_remove(site, "eMMO"),
         genome = as.numeric(genome)) %>%
  select(site, genome, Completeness, Contamination)

genome_data <- read_csv("../data/genome_data.csv") %>%
  mutate(check = if_else(call == "ambiguous" & str_detect(category, "Metal reduction") & str_detect(category, "iron_reduction"), "fegenie", call),
         category = if_else(call == "ambiguous" & str_detect(category, "Metal reduction"), str_split(category, ",")[[1]][1], category),
         call = check) %>% #re-annotate mismatching annotations as "ambiguous
  select(-check) %>%
  group_by(site, genome, call, category, gene_function, Gene.abbreviation, Gene.name) %>%
  summarise(rel_hits = sum(rel_hits), hits = sum(hits)) 

MAG_gene_counts <- read_delim("../data/metabolic/genomes/geneCounts.txt", delim = "\t", col_names = F) %>%
  separate(X1, c("id", "gene_counts"), sep = " ") %>%
  separate(id, c("site", "genome")) %>%
  mutate(site = str_remove(site, "eMMO"))

MAG_annotations <- genome_data %>%
  group_by(call, site) %>%
  summarise(hits = sum(hits)) %>%
  group_by(site) %>%
  mutate(total_hits = sum(hits)) %>%
  left_join(MAG_gene_counts %>% group_by(site) %>% summarise(gene_counts = sum(as.numeric(gene_counts)))) %>%
  group_by(call, site) %>%
  summarise(prop_hits = hits/total_hits,
            rel_hits = hits/as.numeric(gene_counts))

MAG_annotations_plot <- MAG_annotations %>%
  pivot_longer(prop_hits:rel_hits, names_to = "name", values_to = "value") %>%
  mutate(name = str_to_title(str_replace_all(name, "_", " ")),
         site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")))) %>%
  ggplot(aes(site, value, fill = call)) +
  geom_bar(stat="identity") +
  facet_wrap(~name, scales = "free") + 
  coord_flip() +
  ggtitle("B. MAG annotations") +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 18))

Fe_cycler_MAGs <- genome_data %>%
  ungroup() %>%
  filter(category %in% c("iron_oxidation", "iron_reduction", "probable_iron_reduction", "possible_iron_oxidation_and_possible_iron_reduction")) %>%
  select(site, genome) %>%
  distinct() %>%
  mutate(fe_cycler = 1)

MAG_stats_table <- genome_data %>%
  group_by(site, genome) %>%
  summarise(hits = sum(hits)) %>%
  left_join(MAG_gene_counts %>% mutate(genome = as.numeric(genome))) %>% 
  left_join(checkm_stats) %>%
  mutate(`80% Complete` = if_else(Completeness >= 80 & Contamination <= 10, 1, 0)) %>%
  left_join(Fe_cycler_MAGs) %>%
  mutate_if(is.numeric,coalesce,0) %>%
  mutate(complete_fe_cycler = if_else(`80% Complete` == 1 & fe_cycler == 1, 1, 0)) %>%
  left_join(MAG_gene_counts %>% mutate(genome = as.numeric(genome))) %>%
  group_by(site) %>%
  summarise(`80% Complete` = sum(`80% Complete`, na.rm = T), 
            fe_cyclers = sum(fe_cycler),
            complete_fe_cyclers = sum(complete_fe_cycler),
            total_genomes = n(),
            total_genes = sum(as.numeric(gene_counts))) %>%
  # mutate(`80% Complete (%)` = `80% Complete`/total_genomes*100,
  #        `% Fe cyclers` = fe_cyclers/total_genomes*100) %>%
  pivot_longer(-site, names_to = "observation", values_to = "value") %>%
  mutate(observation = str_to_title(str_replace_all(observation, "_", " ")),
         site = factor(site, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC"))) %>%
  ggplot(aes(x=observation,y=site, label = value)) + 
  geom_tile(fill = "white") + geom_text(colour = "black") +
  facet_grid(vars(site),vars(observation), scales = "free", labeller = label_wrap_gen(width = 2, multi_line = TRUE)) +
  ggtitle("D. MAG Stats") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        axis.ticks=element_blank(),
        text = element_text(size = 18))

#write supp table 4
# genome_data %>%
#   group_by(site, genome) %>%
#   summarise(hits = sum(hits)) %>%
#   left_join(MAG_gene_counts %>% mutate(genome = as.numeric(genome))) %>% 
#   left_join(checkm_stats) %>%
#   left_join(taxonomy %>% mutate(genome = as.numeric(genome))) %>%
#   write_csv("../data/supp_table4.csv")

  # plot data ---------------------------------------------------------------

stat_bar_plots <- plot_grid(annotations_plot, MAG_annotations_plot, align = "v", axis = "l")
stat_tables <- plot_grid(metagenome_stats_table, MAG_stats_table, align = "v", axis = "l")
plot_grid(stat_bar_plots, stat_tables, nrow = 2, align = "h")
