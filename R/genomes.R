pacman::p_load(tidyverse, heatmaply, viridis, gplots, dendextend, vegan, gplots, Heatplus, scales, cowplot)

metadata <- read_csv("../data/metadata.csv")
taxonomy <- read_csv("../data/taxonomy.csv") %>%
  mutate(genome = as.character(genome))
genome_data <- read_csv("../data/genome_data.csv") %>%
  mutate(check = if_else(call == "ambiguous" & str_detect(category, "Metal reduction") & str_detect(category, "iron_reduction"), "fegenie", call),
         category = if_else(call == "ambiguous" & str_detect(category, "Metal reduction"), str_split(category, ",")[[1]][1], category),
         call = check) %>% #re-annotate mismatching annotations as "ambiguous
  select(-check) %>%
  group_by(site, genome, call, category, gene_function, Gene.abbreviation, Gene.name) %>%
  summarise(rel_hits = sum(rel_hits), hits = sum(hits))

# categorize metabolic pathways -------------------------------------------

pathway_categories <- read_csv("../data/pathway_categories.csv")
categories <- pathway_categories %>%
  select(category, pathway) %>%
  filter(!is.na(category))

gene_functions <- pathway_categories %>%
  select(gene_function, pathway) %>%
  filter(!is.na(gene_function))

genome_fe_pathways <- genome_data %>%
  ungroup() %>%
  inner_join(gene_functions) %>%
  bind_rows(genome_data %>%
              ungroup() %>%
              inner_join(categories)) %>%
  #mutate(Gene.abbreviation = if_else(is.na(Gene.abbreviation), gene_function, Gene.abbreviation)) %>%
  group_by(site, genome, pathway) %>%
  summarise(rel_hits = sum(rel_hits)) %>%
  left_join(taxonomy) %>%
  mutate(taxa = if_else(!str_detect(genus, "(?=.*[0-9]).*") & !is.na(genus), genus,
                        if_else(!str_detect(family, "(?=.*[0-9]).*") & !is.na(family), family,
                                if_else(!str_detect(order, "(?=.*[0-9]).*") & !is.na(order), order,
                                        if_else(!str_detect(class, "(?=.*[0-9]).*") & !is.na(class), class,
                                                if_else(!is.na(phylum), phylum,
                                                        "Unclassified")))))) %>%
  #unite(col = "taxonomy",  phylum:species, na.rm=TRUE, sep = ";") %>%
  pivot_wider(names_from = pathway, values_from = rel_hits, values_fill = list(rel_hits = 0)) %>%
  group_by(site, genome, phylum, taxa) %>%
  filter_at(vars(iron_oxidation, iron_reduction, probable_iron_reduction, possible_iron_oxidation_and_possible_iron_reduction), any_vars(. > 0)) %>%
  left_join(metadata %>% select(site, genome, Completeness, Contamination)) %>%
  #filter(Completeness >= 80) %>%
  select_if(~ !is.numeric(.) || sum(., na.rm = T) != 0)
  #filter_at(vars(iron_oxidation, iron_reduction, probable_iron_reduction, possible_iron_oxidation_and_possible_iron_reduction), any_vars(. > 0))

pathway_heatmap <- genome_fe_pathways %>%
  select(-domain, -class, -order, -family, -genus, -species, -Completeness, -Contamination) %>%
  pivot_longer(cols = -site:-taxa,  names_to = "pathway", values_to = "rel_hits") %>%
  #filter(pathway != "oxygen reduction") %>%
  mutate(id = paste(site, taxa, sep = "_")) %>%
  #mutate(rel_hits = if_else(rel_hits > 0, 1, 0)) %>%
  ggplot(aes(pathway, reorder(id, rel_hits), fill = rel_hits)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(panel.spacing.y=unit(1,"lines"), 
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_blank()) +
  facet_grid(rows = vars(site), scales = "free_y", space = "free_y")


# clustered heatmap -------------------------------------------------------

genome_pathway_matrix <- genome_fe_pathways %>%
  ungroup() %>%
  mutate(id = paste(site, genome, taxa, sep = "_")) %>%
  #select(-site, -genome, -domain, -taxonomy) %>%
  select(-site, -genome, -domain, -phylum, -class, -order, -family, -genus, -species, -taxa, -`oxygen reduction`, -Completeness, -Contamination) %>%
  #mutate_if(is.numeric, ~ replace(., . > 0, 1)) %>%
  column_to_rownames("id")

#heatmap
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
data.dist <- vegdist(genome_pathway_matrix, method = "bray")
row.clus <- hclust(data.dist, "aver")

data.dist.g <- vegdist(t(genome_pathway_matrix), method = "bray")
col.clus <- hclust(data.dist.g, "aver")

taxa_colors <- genome_fe_pathways %>%
  ungroup() %>%
  select(phylum) %>%
  distinct() %>%
  mutate(color = randomcoloR::distinctColorPalette(16)) 

taxa_color_vec <- data.frame(id = rownames(genome_pathway_matrix)) %>%
  separate(id, into = c("site", "genome", "taxa")) %>%
  left_join(taxonomy %>% select(site, genome, phylum) %>% mutate(genome = as.character(genome))) %>%
  left_join(taxa_colors) %>%
  select(color) %>%
  pull()


heatmap.2(as.matrix(genome_pathway_matrix), 
          Rowv = as.dendrogram(row.clus), 
          Colv = as.dendrogram(col.clus), 
          col = viridis_pal(), 
          margins = c(10, 8),
          RowSideColors = taxa_color_vec,
          trace = "none", 
          density.info = "none",
          keysize = 1,
          key.par=list(mar=c(3,5,4,5))
) # this makes the colour-key legend a little thinner



legend("topright",legend=taxa_colors$phylum, 
       fill=taxa_colors$color, cex=0.4, box.lty=0)



# fe cycling taxonomy -----------------------------------------------------
# taxa_colors <- read_csv("../data/silva138_phylum_colors.csv") %>%
#   mutate(phylum = if_else(phylum == "Cyanobacteria", "Cyanobacteriota", phylum),
#          phylum = if_else(phylum == "Chloroflexi", "Chloroflexota", phylum),
#          phylum = if_else(phylum == "Hadesarchaeaeota", "Hadesarchaeota", phylum),
#          phylum = if_else(phylum == "Omnitrophica", "Omnitrophota", phylum),
#          phylum = if_else(phylum == "AABM512524", "AABM5-125-24", phylum)) %>%
#   arrange(order)
# 
# taxa_color_palette <- taxa_colors$hex.color
# names(taxa_color_palette) <- taxa_colors$phylum
# genome_fe_cycling <- genome_data %>%
#   filter(call == "fegenie") %>%
#   #filter(category %in% c("iron_oxidation", "iron_reduction", "probable_iron_reduction", "possible_iron_oxidation_and_possible_iron_reduction")) %>%
#   mutate(gene_function = if_else(category %in% c("iron_reduction"), 
#                                  str_extract(gene_function, "[^_]+"), gene_function),
#          gene_function = if_else(category %in% c("iron_reduction") & str_detect(gene_function, "(?=.*[0-9]).*"),
#                                  "other", gene_function),
#          gene_function = if_else(!category %in% c("iron_oxidation", "iron_reduction"), category, gene_function)) %>%
#   group_by(site, genome, category, gene_function) %>%
#   summarise(hits = sum(hits), rel_hits = sum(rel_hits)) %>%
#   mutate(gene_function = factor(gene_function, levels = names(fe_cycling_palette)),
#          type = if_else(category %in% c("iron_oxidation", "iron_reduction", "probable_iron_reduction", "possible_iron_oxidation_and_possible_iron_reduction"), "Energy Metabolism", "Housekeeping")) %>%
#   left_join(taxonomy) %>%
#   mutate(taxa = if_else(str_detect(phylum, "Proteobacteria") & !is.na(class), class,
#                                                 if_else(!is.na(phylum), phylum,
#                                                         "Unassigned")),
#          taxa = if_else(str_detect(taxa, "Firmicutes"), "Firmicutes", taxa),
#          taxa = if_else(str_detect(taxa, "Goldbacteria"), "Goldbacteria", taxa),
#          taxa = if_else(str_detect(taxa, "JdFR-18"), "Archaea", taxa),
#          taxa = factor(taxa, levels = taxa_colors$phylum)) %>%
#   group_by(site, taxa, category, gene_function, type) %>%
#   summarise(`Number of Genomes` = n()) %>%
#   left_join(MAG_gene_counts %>% group_by(site) %>% summarise(n_genomes = n())) %>%
#   mutate(site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC"))))

taxa_function <- genome_data %>%
  ungroup() %>%
  filter(category %in% c("iron_oxidation", "iron_reduction")) %>%
  select(site, genome, category) %>%
  distinct() %>%
  left_join(taxonomy) %>%
  mutate(taxa = if_else(!str_detect(family, "(?=.*[0-9]).*") & !is.na(family), family,
                        if_else(!str_detect(order, "(?=.*[0-9]).*") & !is.na(order), order,
                                if_else(!str_detect(class, "(?=.*[0-9]).*") & !is.na(class), class,
                                        if_else(!is.na(phylum), phylum,
                                                "Unclassified")))),
         taxa = if_else(class %in% c("Thermodesulfovibrionia", "Desulfobulbia", "Nitrospiria"), class, taxa),
         taxa = if_else(order %in% c("Ignavibacteriales", "Desulfobacterales", "Obscuribacterales"), order, taxa)) %>%
  mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum)) %>%
  group_by(site, phylum, taxa, category) %>%
  summarise(`Number of Genomes` = n()) %>%
  left_join(MAG_gene_counts %>% group_by(site) %>% summarise(n_genomes = n())) %>%
  mutate(`% of Genomes` = `Number of Genomes`/n_genomes*100) 

taxa_genes <- genome_data %>%
  ungroup() %>%
  filter(category %in% c("iron_oxidation", "iron_reduction")) %>%
  select(site, genome, category, gene_function) %>%
  distinct() %>%
  left_join(taxonomy %>% mutate(genome = as.numeric(genome))) %>%
  mutate(taxa = if_else(!str_detect(family, "(?=.*[0-9]).*") & !is.na(family), family,
                                if_else(!str_detect(order, "(?=.*[0-9]).*") & !is.na(order), order,
                                        if_else(!str_detect(class, "(?=.*[0-9]).*") & !is.na(class), class,
                                                if_else(!is.na(phylum), phylum,
                                                        "Unclassified")))),
         taxa = if_else(class %in% c("Thermodesulfovibrionia", "Desulfobulbia", "Nitrospiria"), class, taxa),
         taxa = if_else(order %in% c("Ignavibacteriales", "Desulfobacterales", "Obscuribacterales"), order, taxa)) %>%
  mutate(gene_function = if_else(category %in% c("iron_reduction"), 
                                 str_extract(gene_function, "[^_]+"), gene_function),
         gene_function = if_else(category %in% c("iron_reduction") & str_detect(gene_function, "(?=.*[0-9]).*"),
                                 "other", gene_function),
         gene_function = if_else(!category %in% c("iron_oxidation", "iron_reduction"), category, gene_function)) %>%
  group_by(site, taxa, category, gene_function) %>%
  summarise(`Number of Genomes` = n()) %>%
  left_join(MAG_gene_counts %>% group_by(site) %>% summarise(n_genomes = n())) %>%
  mutate(`% of Genomes` = `Number of Genomes`/n_genomes*100) 

fe_cyclers <- taxa_function %>% 
  ungroup() %>%  
  group_by(taxa) %>%
  summarise(`% of Genomes` = sum(`% of Genomes`)) %>%
  arrange(desc(`% of Genomes`)) %>%
  select(taxa) %>% distinct() %>% pull()

taxa_rainbow_palette <- colorRampPalette(c("pink","red", "orange", "yellow", "springgreen","royalblue", "purple", "gray", "#2e2d2d"))(length(fe_cyclers))
names(taxa_rainbow_palette) <- c(fe_cyclers[2:length(fe_cyclers)], fe_cyclers[1])
  

taxa_function_plot <- taxa_function %>%
  mutate(site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC"))),
         taxa = factor(taxa, levels = names(taxa_rainbow_palette))) %>%
  ggplot(aes(site, `% of Genomes`, fill = taxa)) +
  geom_bar(stat ="identity") +
  coord_flip() +
  scale_fill_manual(values = taxa_rainbow_palette, name = "Taxa") +
  theme(axis.title.y = element_blank()) +
  facet_grid(rows = vars(category),labeller = label_wrap_gen(multi_line = TRUE)) +
  theme(text = element_text(size = 18),
        strip.text = element_text(face="bold")) +
  ggtitle("A. Iron Cycling Taxa")

taxa_gene_feox <- taxa_genes %>%
  filter(category == "iron_oxidation") %>%
  mutate(gene_function = str_replace_all(gene_function, "repCluster", "Cluster"),
         gene_function = str_replace_all(gene_function, "_", " "),
         site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC"))),
         taxa = factor(taxa, levels = names(taxa_rainbow_palette))) %>%
  ggplot(aes(site, `% of Genomes`, fill = taxa)) +
  geom_bar(stat ="identity") +
  coord_flip() +
  scale_fill_manual(values = taxa_rainbow_palette, name = "Taxa") +
  theme(axis.title.y = element_blank()) +
  facet_grid(rows = vars(gene_function), labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
  theme(text = element_text(size = 18),
        strip.text = element_text(face="bold"),
        legend.position = "none")  +
  ggtitle("B. Iron Oxidation Genes")

taxa_gene_fered <- taxa_genes %>%
  filter(category == "iron_reduction") %>%
  mutate(gene_function = str_replace_all(gene_function, "repCluster", "Cluster"),
         gene_function = str_replace_all(gene_function, "_", " "),
         site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC"))),
         taxa = factor(taxa, levels = names(taxa_rainbow_palette))) %>%
  ggplot(aes(site, `% of Genomes`, fill = taxa)) +
  geom_bar(stat ="identity") +
  coord_flip() +
  scale_fill_manual(values = taxa_rainbow_palette, name = "Taxa") +
  theme(axis.title.y = element_blank()) +
  facet_grid(rows = vars(gene_function), labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
  theme(text = element_text(size = 18),
        strip.text = element_text(face="bold"),
        legend.position = "none")  +
  ggtitle("B. Iron Reduction Genes")

plot_grid(taxa_function_plot, 
          plot_grid(taxa_gene_feox, taxa_gene_fered, align = "v", axis = "l", rel_widths = c(2, 2), ncol = 2), nrow = 2, rel_heights = c(1, 2))


# genome_function_plot_housekeeping <- genome_fe_cycling %>%
#   filter(type == "Housekeeping") %>%
#   group_by(site, taxa, category) %>%
#   summarise(`Number of Genomes` = sum(`Number of Genomes`)) %>%
#   ggplot(aes(site, `Number of Genomes`, fill = taxa)) +
#   geom_bar(stat ="identity", position = "fill") +
#   coord_flip() +
#   scale_fill_manual(values = taxa_color_palette, name = "Taxa") +
#   theme(axis.title.y = element_blank()) +
#   facet_grid(rows = vars(category),labeller = label_wrap_gen(multi_line = TRUE)) +
#   theme(strip.text = element_text(face="bold", size=6)) +
#   ggtitle("B. Housekeeping") +
#   guides(fill = guide_legend(ncol = 1))
# 
# plot_grid(genome_function_plot, genome_function_plot_housekeeping, align = "v", axis = "l", rel_widths = c(2, 2.5))
# 

# genome_energy_plot <- genome_fe_cycling %>%
#   filter(type == "Energy Metabolism") %>%
#   ggplot(aes(site, `Number of Genomes`, fill = taxa)) +
#   geom_bar(stat ="identity", position = "fill") +
#   coord_flip() +
#   scale_fill_manual(values = taxa_color_palette, name = "Taxa") +
#   theme(axis.title.y = element_blank()) +
#   facet_grid(rows = vars(gene_function), labeller = label_wrap_gen(multi_line = TRUE)) +
#   theme(legend.position = "none",
#         strip.text = element_text(face="bold", size=6)) +
#   ggtitle("A. Energy Metabolism Genes")
# 
# genome_housekeeping_plot <- genome_fe_cycling %>%
#   filter(type == "Housekeeping") %>%
#   ggplot(aes(site, `Number of Genomes`, fill = taxa)) +
#   geom_bar(stat ="identity", position = "fill") +
#   coord_flip() +
#   scale_fill_manual(values = taxa_color_palette, name = "Taxa") +
#   theme(axis.title.y = element_blank(),
#         strip.text = element_text(face="bold", size=6)) +
#   facet_grid(rows = vars(gene_function), labeller = label_wrap_gen(multi_line = TRUE)) +
#   ggtitle("B. Housekeeping Genes") +
#   guides(fill = guide_legend(ncol = 1))
# 
# plot_grid(genome_energy_plot, genome_housekeeping_plot, align = "v", axis = "l", rel_widths = c(2, 2.5))



# bubble plot -------------------------------------------------------------

genome_bubble_plot <- genome_data %>%
  filter(call == "fegenie") %>%
  left_join(taxonomy %>% mutate(genome = as.numeric(genome))) %>%
  left_join(metadata %>% select(site, genome, Completeness)) %>%
  mutate(taxa = if_else(!str_detect(family, "(?=.*[0-9]).*") & !is.na(family), family,
                        if_else(!str_detect(order, "(?=.*[0-9]).*") & !is.na(order), order,
                                if_else(!str_detect(class, "(?=.*[0-9]).*") & !is.na(class), class,
                                        if_else(!is.na(phylum), phylum,
                                                "Unclassified")))),
         taxa = if_else(class %in% c("Thermodesulfovibrionia", "Desulfobulbia", "Nitrospiria"), class, taxa),
         taxa = if_else(order %in% c("Ignavibacteriales", "Desulfobacterales", "Obscuribacterales"), order, taxa)) %>%
  mutate(id = paste(site, genome, taxa, sep = "_")) %>%
  mutate(category = if_else(str_detect(category, "aquisition"), "iron acquisition", category),
         category = if_else(str_detect(category, "storage"), "iron storage", category)) %>%
  group_by(site, id, taxa, Completeness, category) %>%
  summarise(rel_hits = sum(rel_hits)) %>%
  mutate(rel_hits = 1) %>%
  pivot_wider(names_from = category, values_from = rel_hits, values_fill = list(rel_hits = 0)) %>%
  mutate(housekeeping = sum(`iron acquisition`, `iron storage`, iron_gene_regulation)) %>%
  group_by(`iron acquisition`, `iron storage`, iron_gene_regulation) %>%
  #summarise(total = n()/515*100)
  ggplot(aes(category, reorder(id, rel_hits), fill = log10(rel_hits), label=taxa)) +
  geom_tile() +
  #scale_size_continuous(breaks = c(2e-05, 5e-05, 1e-04, 5e-04, 1e-03), name = "% metagenome") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  guides(col = guide_legend(ncol = 1)) +
  scale_fill_viridis_c() +
  #guides(fill=guide_legend(title="Category")) +
  theme(axis.title.x=ggplot2::element_blank(), 
        axis.title.y=ggplot2::element_blank(),
        #legend.position = "none",
        #strip.background = element_blank(), 
        panel.spacing = unit(0,"line"), 
        panel.border = element_rect(size = 0.25, color = "black"),
        strip.text.y = element_text(angle = 180, size=8, lineheight=1)) +
  facet_grid(rows = vars(site), scales = "free")

plotly::ggplotly(genome_bubble_plot)
