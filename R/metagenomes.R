pacman::p_load(tidyverse, heatmaply, viridis, gplots, dendextend, vegan, gplots, Heatplus, scales)

metagenome_data <- read_csv(paste0(write_path, "metagenome_data.csv")) %>%
  mutate(check = if_else(call == "ambiguous" & str_detect(category, "Metal reduction") & str_detect(category, "iron_reduction"), "fegenie", call),
         category = if_else(call == "ambiguous" & str_detect(category, "Metal reduction"), str_split(category, ",")[[1]][1], category),
         call = check) %>% #re-annotate mismatching annotations as "ambiguous
  select(-check) %>%
  group_by(site, genome, call, category, gene_function, Gene.abbreviation, Gene.name) %>%
  summarise(rel_hits = sum(rel_hits), hits = sum(hits)) %>%
  mutate(date = "Apr 2018")

momper_2017_metagenome_data <- read_csv(paste0(data_path, "momper2017_fegenie_metagenome_data.csv"))

# categorize metabolic pathways -------------------------------------------

pathway_categories <- read_csv(paste0(data_path, "pathway_categories.csv")) %>%
  filter(!gene_function %in% c("Nitrite reduction", "Nitrite reduction to ammonia", "Nitrite oxidation")) %>%
  filter(!category %in% c("possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction"))
categories <- pathway_categories %>%
  select(category, pathway) %>%
  filter(!is.na(category))

gene_functions <- pathway_categories %>%
  select(gene_function, pathway) %>%
  filter(!is.na(gene_function))

fe_pathways <- metagenome_data %>%
  ungroup() %>%
  inner_join(gene_functions) %>%
  bind_rows(metagenome_data %>%
              ungroup() %>%
              inner_join(categories)) %>%
  #mutate(Gene.abbreviation = if_else(is.na(Gene.abbreviation), gene_function, Gene.abbreviation)) %>%
  group_by(site, pathway) %>%
  summarise(rel_hits = sum(rel_hits))

#write supp table 1
metagenome_data %>%
 ungroup() %>%
 inner_join(gene_functions) %>%
 bind_rows(metagenome_data %>%
             ungroup() %>%
             inner_join(categories)) %>% select(pathway, call, category, gene_function, Gene.abbreviation) %>% distinct() %>% arrange(pathway, call, category, gene_function) %>% write_csv(paste0(write_path, "supp_table1.csv"))

fe_pathways_plot <- fe_pathways %>%
  ggplot(aes(rel_hits, site, color = pathway, group = pathway)) +
  geom_point() +
  scale_x_log10() +
  facet_grid(cols = vars(pathway))

metagenome_pathway_heatmap <- fe_pathways %>%
  pivot_wider(names_from = pathway, values_from = rel_hits, values_fill = list(rel_hits= 0)) %>%
  select_if(~ !is.numeric(.) || sum(., na.rm = T) != 0) %>%
  pivot_longer(-site, names_to = "pathway", values_to = "rel_hits") %>%
  #filter(!site %in% c("WC")) %>%
  #filter(!pathway %in% c("oxygen reduction")) %>%
  #filter(!pathway %in% c("oxygen reduction", "Carbon fixation")) %>%
  mutate(pathway = str_to_sentence(str_replace_all(pathway, "_", " "))) %>%
  mutate(site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")))) %>%
  #mutate(rel_hits = if_else(rel_hits > 0, 1, 0)) %>%
  ggplot(aes(reorder(pathway, -rel_hits), site, fill = log10(rel_hits))) +
  geom_tile() +
  scale_fill_viridis() +
  labs(fill='Log % Metagenome') +
  theme(panel.spacing.y=unit(1,"lines"), 
        #axis.text.y = element_text(size = 4),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_blank(),
        text = element_text(size = 18)) 

metagenome_all_pathways <- metagenome_data %>%
  filter(call != "ambiguous") %>%
  filter(!category %in% c("Sulfur cycling enzymes (detailed)", "iron_aquisition-heme_oxygenase",
                          "iron_aquisition-heme_transport", "iron_aquisition-iron_transport", "iron_storage",
                          "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport_potential",
                          "iron_gene_regulation", "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction")) %>%
  mutate(category = if_else(category %in% c("iron_oxidation", "iron_reduction", "Metal reduction"), "Iron cycling", category)) %>%
  #filter(category %in% c("Sulfur cycling", "Nitrogen cycling", "Iron cycling", "As cycling", "Methane metabolism", "Hydrogenases")) %>%
  group_by(site, category) %>%
  summarise(rel_hits = sum(rel_hits))

metagenome_all_pathways_heatmap <- metagenome_all_pathways %>%
  mutate(category = str_to_sentence(str_replace_all(category, "_", " "))) %>%
  #filter(!site %in% c("WC", "SW")) %>%
  pivot_wider(names_from = category, values_from = rel_hits, values_fill = list(rel_hits= 0)) %>%
  select_if(~ !is.numeric(.) || sum(., na.rm = T) != 0) %>%
  pivot_longer(-site, names_to = "pathway", values_to = "rel_hits") %>%
  mutate(site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")))) %>%
  #mutate(rel_hits = if_else(rel_hits > 0, 1, 0)) %>%
  ggplot(aes(reorder(pathway, -rel_hits), site, fill = log10(rel_hits))) +
  geom_tile() +
  scale_fill_viridis() +
  labs(fill='Log % Metagenome') +
  theme(panel.spacing.y=unit(1,"lines"), 
        #axis.text.y = element_text(size = 4),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_blank(),
        text = element_text(size = 18)) 

# iron gene diversity -----------------------------------------------------

fe_gene_diversity <- metagenome_data %>%
  filter(call == "fegenie") %>%
  group_by(site, category) %>%
  summarise(rel_hits = sum(rel_hits, na.rm = T), n_genes = n()) 

diversity_plot <- fe_gene_diversity %>%
  ggplot(aes(rel_hits, n_genes, color = category, shape = site, group = category)) +
  geom_point() +
  scale_x_log10() +
  scale_shape_manual(values = c(0,1,2,15,16,17, 4, 8)) +
  facet_wrap(~category, scales = "free_y")


# heat map clustering on iron genes ---------------------------------------

summary_mat <- metagenome_data %>%
  ungroup() %>%
  filter(call == "fegenie") %>%
  mutate(category = recode(category, igr = "iron_gene_regulation",
                           iron_gene_regulation = "igr",
                           `iron_aquisition-heme_transport` = "iaht",
                           `iron_aquisition-siderophore_transport_potential` = "iastp",
                           iron_storage = "is",
                           `iron_aquisition-iron_transport` = "iait",
                           iron_reduction = "ir", 
                           iron_oxidation = "io", 
                           `iron_aquisition-siderophore_synthesis` = "iass",
                           possible_iron_oxidation_and_possible_iron_reduction = "popr",
                           `iron_aquisition-heme_oxygenase` = "iaho",
                           probable_iron_reduction = "pr")) %>%
  filter(category %in% c("ir", "io", "popr", "pr"))

gene_ids <- summary_mat %>%
  ungroup() %>%
  select(category, gene_function) %>%
  distinct() %>%
  arrange(category, gene_function) %>%
  group_by(category) %>%
  mutate(id = 1:n(),
         gene_id= paste(category, id, sep = "_"))



#write data to Supplementary table
gene_ids %>% 
  mutate(category = recode(category, io = "iron oxidation",
                           ir = "iron reduction",
                           popr = "probable iron oxidation or reduction",
                           pr = "probable iron reduction")) %>%
  select(gene_id, category, gene_function) %>% write_csv(paste0(write_path, "supp_table4.csv"))

summary_matrix <- summary_mat %>%
  left_join(gene_ids %>% select(category, gene_function, gene_id)) %>%
  select(site, gene_id, rel_hits) %>%
  pivot_wider(id_cols = site, names_from = gene_id, values_from = rel_hits, values_fill = list(rel_hits = 0)) %>%
  add_row(site = "WC") %>%
  mutate_if(is.numeric,coalesce,0) %>%
  column_to_rownames("site") 

#heatmap
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
data.dist <- vegdist(summary_matrix, method = "bray")
row.clus <- hclust(data.dist, "aver")

data.dist.g <- vegdist(t(summary_matrix), method = "bray")
col.clus <- hclust(data.dist.g, "aver")

cat_colors <- metagenome_data %>%
  ungroup() %>%
  filter(call == "fegenie") %>% 
  select(category) %>%
  distinct() %>%
  mutate(color = hue_pal()(11)) %>%
  mutate(abbr = recode(category, igr = "iron_gene_regulation",
                           iron_gene_regulation = "igr",
                           `iron_aquisition-heme_transport` = "iaht",
                           `iron_aquisition-siderophore_transport_potential` = "iastp",
                           iron_storage = "is",
                           `iron_aquisition-iron_transport` = "iait",
                           iron_reduction = "ir", 
                           iron_oxidation = "io", 
                           `iron_aquisition-siderophore_synthesis` = "iass",
                           possible_iron_oxidation_and_possible_iron_reduction = "popr",
                           `iron_aquisition-heme_oxygenase` = "iaho",
                           probable_iron_reduction = "pr"))

cat_color_vec <- data.frame(gene_id = colnames(summary_matrix)) %>%
  mutate(id = gene_id) %>%
  separate(gene_id, into = c("abbr", "gene_id")) %>%
  left_join(cat_colors) %>%
  select(color) %>%
  pull()



# iron cycling breakdown --------------------------------------------------
fe_cycling_palette <- c(colorRampPalette(c("springgreen","royalblue", "purple"))(6),
                        colorRampPalette(c("brown", "red","orange", "yellow"))(8),
                        colorRampPalette(c("gray", "black"))(2),
                        colorRampPalette(c("#faa7e8", "#d6faa7", "#a7faf2", "#cabbfc"))(7))
names(fe_cycling_palette) <- str_replace_all(c("Cyc1","Cyc2_repCluster1","Cyc2_repCluster2","Cyc2_repCluster3","FoxY","FoxE",
                               "DFE","OmcF","OmcS","OmcZ",
                               "MtrA","MtrB","MtrC","other",
                               "possible_iron_oxidation_and_possible_iron_reduction","probable_iron_reduction",
                               "iron_aquisition-heme_oxygenase","iron_aquisition-heme_transport","iron_aquisition-iron_transport","iron_aquisition-siderophore_synthesis",              
                               "iron_aquisition-siderophore_transport_potential","iron_gene_regulation","iron_storage"), "_", " ")

metagenome_fe_cycling <- metagenome_data %>%
  filter(call == "fegenie") %>%
  bind_rows(momper_2017_metagenome_data) %>%
  #filter(category %in% c("iron_oxidation", "iron_reduction", "probable_iron_reduction", "possible_iron_oxidation_and_possible_iron_reduction")) %>%
  mutate(gene_function = if_else(category %in% c("iron_reduction"), 
                                 str_extract(gene_function, "[^_]+"), gene_function),
         gene_function = if_else(category %in% c("iron_reduction") & str_detect(gene_function, "(?=.*[0-9]).*"),
                                 "other", gene_function),
         gene_function = if_else(!category %in% c("iron_oxidation", "iron_reduction"), category, gene_function),
         gene_function = str_replace_all(gene_function, "_", " ")) %>%
  group_by(site, date, category, gene_function) %>%
  summarise(hits = sum(hits), rel_hits = sum(rel_hits)) %>%
  mutate(gene_function = factor(gene_function, levels = names(fe_cycling_palette)),
         site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC", "DuselD"))),
         type = if_else(category %in% c("iron_oxidation", "iron_reduction", "probable_iron_reduction", "possible_iron_oxidation_and_possible_iron_reduction"), "Energy Metabolism", "Housekeeping")) %>%
  ggplot(aes(site, rel_hits, fill = gene_function)) +
  geom_bar(stat ="identity") +
  coord_flip() +
  scale_fill_manual(values = fe_cycling_palette, name = "Gene") +
  guides(fill = guide_legend(ncol = 1)) +
  ylab("% Metagenome") +
  theme(axis.title.y = element_blank()) +
  facet_grid(cols = vars(type), rows = vars(date), scales = "free", space = "free_y") +
  theme(text = element_text(size = 18)) 
  #facet_wrap(~type, scales = "free_x")



#housekeeping breakdown
housekeeping_genes <- metagenome_data %>%
  filter(call == "fegenie") %>%
  filter(!category %in% c("iron_oxidation", "iron_reduction", "probable_iron_reduction", "possible_iron_oxidation_and_possible_iron_reduction")) %>%
  arrange(desc(rel_hits)) %>%
  ungroup() %>%
  select(category, gene_function) %>%
  distinct()
  


