pacman::p_load(tidyverse, readxl)

jgi_data <- read_delim("../data/JGI_April2018.tsv", delim = "\t") %>%
  filter(product == "OmcA/MtrC family decaheme c-type cytochrome") %>%
  select(where(~!all(is.na(.)))) %>%
  pivot_longer(-product:-source, names_to = "IMG_AP", values_to = "hits")

jgi_annotations <- jgi_data %>%
  select(product) %>%
  distinct()

jgi_metadata_list <- "../data/JGI_April2018_metadata.xlsx" %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = path)

JGI_metadata <- reduce(jgi_metadata_list, full_join)

jgi_hits <- jgi_data %>%
  left_join(JGI_metadata)
