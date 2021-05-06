pacman::p_load(tidyverse, janitor)


momper2017_gene_count <- read_delim(paste0(data_path, "fegenie/metagenomes/momper2017/geneCounts.txt"), delim = "\t", col_names = F) %>%
  separate(X1, c("id", "gene_counts"), sep = "  ") %>%
  separate(id, c("site", "genome"))  %>%
  mutate(gene_counts = as.numeric(gene_counts))

momper2017_mag_metadata <- read_csv(paste0(data_path, "fegenie/genomes/momper2017/MAG_metadata.csv")) %>%
  remove_empty()


read_fegenie <- function(type){
  files <- list.files(paste0(paste0(data_path, "fegenie/"), type, "/momper2017"), full.names = T, pattern = ".*geneSummary.csv")
  read_file <- function(file){
    site <- str_extract(file, "(?<=momper2017/)(.*)(?=_FeGenie)")
    file %>%
    read_csv(skip = 1) %>%
    mutate(site = site,
           id = str_remove(`genome/assembly`, "[.]fa"),
           fegenie = 1,
           date = "Oct 2013") %>%
    rename(gene_function = "HMM") %>%
    select(site, id, date, orf, fegenie, category, gene_function)
  }
  
  fegenie_list = lapply(files, read_file)
  fegenie_data <- reduce(fegenie_list, full_join)
}

#read fegenie data
momper2017_meta_fegenie <- read_fegenie("metagenomes") %>%
  group_by(site, date, category, gene_function) %>%
  summarise(hits = n()) %>%
  left_join(momper2017_gene_count) %>%
  mutate(rel_hits = hits/gene_counts*100)

write_csv(momper2017_meta_fegenie, paste0(write_path, "momper2017_fegenie_metagenome_data.csv"))

momper2017_genomes_fegenie <- read_fegenie("genomes") %>%
  left_join(momper2017_mag_metadata) 


# plot MAG annotations ----------------------------------------------------

momper2017_taxa_function <- momper2017_genomes_fegenie %>%
  filter(category %in% c("iron_oxidation", "iron_reduction")) %>%
  select(id, Taxon, category) %>%
  distinct() %>%
  separate(Taxon, c("phylum", "taxa"), "; ") %>%
  mutate(taxa = if_else(is.na(taxa), phylum, taxa),
         taxa = if_else(taxa %in% c("OP3", "Omnitrophica"), "Omnitrophota", taxa),
         taxa = if_else(taxa %in% c("Desulfurivibrio", "Desulfobacteraceae"), "Desulfobacterales", taxa),
         taxa = if_else(taxa %in% c("Myxococcales"), "Myxococcota", taxa),
         taxa = if_else(taxa %in% c("Nitrospiraceae"), "Nitrospiria", taxa)) %>%
  group_by(taxa, category) %>%
  summarise(`Number of Genomes` = n()) %>%
  mutate(`% of Genomes` = `Number of Genomes`/74*100,
         site = "D6/DuselD Oct 2013") 

momper2017_fe_cyclers_all <- momper2017_taxa_function %>% ungroup() %>%  select(taxa) %>% distinct() %>% pull()
momper2017_fe_cyclers <- momper2017_fe_cyclers_all[which(!momper2017_fe_cyclers_all %in% fe_cyclers)]

taxa_rainbow_palette <- c(colorRampPalette(c("pink","red", "orange", "yellow", "springgreen","royalblue", "purple", "gray", "#2e2d2d"))(length(fe_cyclers)),
                          colorRampPalette(c("#852222", "#853e0c", "#a19e1f", "#14610d", "#0f305c"))(length(momper2017_fe_cyclers)))
names(taxa_rainbow_palette) <- c(fe_cyclers[2:length(fe_cyclers)], fe_cyclers[1], momper2017_fe_cyclers)

momper2017_taxa_genes <- momper2017_genomes_fegenie %>%
  filter(category %in% c("iron_oxidation", "iron_reduction")) %>%
  select(id, Taxon, category, gene_function) %>%
  distinct() %>%
  separate(Taxon, c("phylum", "taxa"), "; ") %>%
  mutate(taxa = if_else(is.na(taxa), phylum, taxa),
         taxa = if_else(taxa %in% c("OP3", "Omnitrophica"), "Omnitrophota", taxa),
         taxa = if_else(taxa %in% c("Desulfurivibrio", "Desulfobacteraceae"), "Desulfobacterales", taxa),
         taxa = if_else(taxa %in% c("Myxococcales"), "Myxococcota", taxa),
         taxa = if_else(taxa %in% c("Nitrospiraceae"), "Nitrospiria", taxa)) %>%
  mutate(gene_function = if_else(category %in% c("iron_reduction"), 
                                 str_extract(gene_function, "[^_]+"), gene_function),
         gene_function = if_else(category %in% c("iron_reduction") & str_detect(gene_function, "(?=.*[0-9]).*"),
                                 "other", gene_function),
         gene_function = if_else(!category %in% c("iron_oxidation", "iron_reduction"), category, gene_function)) %>%
  group_by(taxa, category, gene_function) %>%
  summarise(`Number of Genomes` = n()) %>%
  mutate(`% of Genomes` = `Number of Genomes`/74*100,
         site = "D6/DuselD Oct 2013") 


momper2017_taxa_function_plot <- momper2017_taxa_function %>%
  mutate(category = str_to_title(str_replace_all(category, "_", " ")),
         taxa = factor(taxa, levels = names(taxa_rainbow_palette))) %>%
  ggplot(aes(site, `% of Genomes`, fill = taxa)) +
  geom_bar(stat ="identity") +
  coord_flip() +
  scale_fill_manual(values = taxa_rainbow_palette, name = "Taxa") +
  theme(axis.title.y = element_blank()) +
  facet_grid(rows = vars(category),labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
  theme(strip.text = element_text(face="bold"),
        text = element_text(size = 18)) +
  ggtitle("A. Iron Cycling Taxa") +
  guides(fill = guide_legend(ncol = 2))

momper2017_taxa_gene_feox <- momper2017_taxa_genes %>%
  filter(category == "iron_oxidation") %>%
  mutate(gene_function = str_replace_all(gene_function, "repCluster", "Cluster"),
         gene_function = str_replace_all(gene_function, "_", " "),
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

momper2017_taxa_gene_fered <- momper2017_taxa_genes %>%
  filter(category == "iron_reduction") %>%
  mutate(gene_function = str_replace_all(gene_function, "repCluster", "Cluster"),
         gene_function = str_replace_all(gene_function, "_", " "),
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


