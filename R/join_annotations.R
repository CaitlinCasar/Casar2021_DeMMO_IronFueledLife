pacman::p_load(tidyverse, readxl, splitstackshape)

data_types <- c("genomes", "metagenomes")

join_annotations <- function(data_type){
  
  path <- paste0(data_path, "fegenie/", data_type)
  files <- list.files(path, full.names = T, pattern = ".*geneSummary.csv")
  
  read_fegenie <- function(file){
    file %>%
    read_csv(skip = 1) %>%
      mutate(id = str_remove(`genome/assembly`, "[.]faa"),
             fegenie = 1) %>%
      select(id, orf, fegenie, category, HMM)
  }
  
  #read fegenie data
  fegenie_list = lapply(files, read_fegenie)
  fegenie_data <- reduce(fegenie_list, full_join)
  
  
  path <- paste0(data_path, "metabolic/", data_type)
  files <- list.files(path, full.names = T, pattern = ".*METABOLIC_result.xlsx")
  
  gene_counts <- read_delim(paste0(path, "/geneCounts.txt"), delim = "\t", col_names = F) %>%
    separate(X1, c("id", "gene_counts"), sep = " ") %>%
    separate(id, c("site", "genome")) %>%
    mutate(site = str_remove(site, "eMMO"))
  
  #read metabolic data
  read_metabolic <- function(file){
    file %>%
    read_xlsx(col_types = "text") %>%
    pivot_longer(-Category:-Hmm.detecting.threshold, names_to = "name", values_to = "value") %>%
    separate(name, c("id", "name"), "[.]") %>%
    pivot_wider(id_cols = Category:id, names_from = name, values_from = value) %>%
    filter(Hit != "0") %>%
    mutate(Hits = str_remove_all(Hits, "None,|,None"),
           Hits = str_remove_all(Hits, paste0(id, "_"))) %>%
    cSplit('Hits', ',') %>%
    pivot_longer(-Category:-Hit, names_to = "hit_id", values_to ="orf", values_drop_na = T) 
  }
  metabolic_list = lapply(files, read_metabolic)
  metabolic_data <- reduce(metabolic_list, full_join)
    
  metabolic_data %>%
    mutate(Hit = 1) %>%
    full_join(fegenie_data) %>%
    mutate(Hit = if_else(is.na(Hit), 0, Hit),
           fegenie = if_else(is.na(fegenie), 0, fegenie),
           call = if_else(Hit > 0, "metabolic", "fegenie"),
           call = if_else(Hit > 0 & fegenie == 1, "ambiguous", call)) %>%
    filter(id != "total") %>%
    separate(id, c("site", "genome")) %>%
    mutate(category = if_else(is.na(category), Category, 
                              if_else(call == "ambiguous", paste(category, Category, sep = ","), category)),
           gene_function = if_else(call == "fegenie", HMM, 
                                   if_else(call == "ambiguous", paste(HMM, Function, sep = ","), Function)),
           site = str_remove(site, "eMMO")) %>%
    group_by(site, genome, call, category, gene_function, Gene.abbreviation, Gene.name) %>%
    summarise(hits = n()) %>%
    left_join(gene_counts) %>%
    mutate(gene_counts = as.numeric(gene_counts),
           rel_hits = hits/gene_counts*100)
}

genome_data <- join_annotations(data_types[1]) %>%
  write_csv(paste0(write_path, "genome_data.csv"))
metagenome_data <- join_annotations(data_types[2])%>%
  write_csv(paste0(write_path, "metagenome_data.csv"))


  
         