pacman::p_load(tidyverse, heatmaply, viridis, gplots, dendextend, vegan, gplots, Heatplus, scales, cowplot)

metadata <- read_csv(paste0(data_path, "metadata.csv")) %>%
  mutate(genome = as.character(genome)) %>%
  select(-Phylum)
taxonomy <- read_csv(paste0(data_path, "taxonomy.csv")) %>%
  mutate(genome = as.character(genome))
genome_data <- read_csv(paste0(write_path, "genome_data.csv")) %>%
  mutate(check = if_else(call == "ambiguous" & str_detect(category, "Metal reduction") & str_detect(category, "iron_reduction"), "fegenie", call),
         category = if_else(call == "ambiguous" & str_detect(category, "Metal reduction"), str_split(category, ",")[[1]][1], category),
         call = check,
         genome = as.character(genome)) %>% #re-annotate mismatching annotations as "ambiguous
  select(-check) %>%
  group_by(site, genome, call, category, gene_function, Gene.abbreviation, Gene.name) %>%
  summarise(rel_hits = sum(rel_hits), hits = sum(hits))


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
  left_join(metadata %>% group_by(site) %>% summarise(n_genomes = n())) %>%
  mutate(`% of Genomes` = `Number of Genomes`/n_genomes*100) 

taxa_genes <- genome_data %>%
  ungroup() %>%
  filter(category %in% c("iron_oxidation", "iron_reduction")) %>%
  select(site, genome, category, gene_function) %>%
  distinct() %>%
  left_join(taxonomy) %>%
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
  left_join(metadata %>% group_by(site) %>% summarise(n_genomes = n())) %>%
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


