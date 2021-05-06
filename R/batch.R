require(pacman)

# set file paths ---------------------------------------------------------------
write_path <- "../../script_test/"
data_path <- "../data/"
write_figures <- paste0(write_path, "figures/")


# join FeGenie and METABOLIC datasets --------------------------------------

source("join_annotations.R")


# calculate data summary stats --------------------------------------------

source("data_stats.R")

## write data stats figure 
pdf(paste0(write_figures, "Supp_Figure1.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
data_stats_plot # print it
dev.off()


# plot metagenome data ----------------------------------------------------

source("metagenomes.R")

## write metagenome-level iron cycling bar plot figure 
pdf(paste0(write_figures, "Figure1.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
metagenome_fe_cycling # print it
dev.off()

## write FeGenie metagenome heatmap figure 
pdf(paste0(write_figures, "Figure2.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
heatmap.2(as.matrix(summary_matrix), 
          Rowv = as.dendrogram(row.clus), 
          Colv = as.dendrogram(col.clus), 
          col = viridis_pal(), 
          margins = c(10, 8),
          ColSideColors = cat_color_vec,
          trace = "none", 
          density.info = "none",
          keysize = 1,
          key.par=list(mar=c(3,5,4,5))
) # this makes the colour-key legend a little thinner

legend("topright",legend=cat_colors$category, 
       fill=cat_colors$color, cex=0.4, box.lty=0)
dev.off()

## write metagenome metabolic pathway heatmap figure 
pdf(paste0(write_figures, "Figure3.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
metagenome_pathway_heatmap
dev.off()

## write metagenome nutrient cycle heatmap figure 
pdf(paste0(write_figures, "Supp_Figure2.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
metagenome_all_pathways_heatmap
dev.off()


# plot MAG data -----------------------------------------------------------

source("genomes.R")

## write MAG-level iron cycling bar plot figure 
pdf(paste0(write_figures, "Figure4.pdf"), height = 11, width = 8.5, paper = "letter")
plot_grid(taxa_function_plot, 
          plot_grid(taxa_gene_feox, taxa_gene_fered, align = "v", axis = "l", rel_widths = c(2, 2), ncol = 2), nrow = 2, rel_heights = c(1, 2))
dev.off()

# geochem data prep for input to thermo models -------------------------------------------------------
source("geochem.R")

#optional plot geochem data
geochem_plot


# plot thermodynamic models -----------------------------------------------
source("thermo.R")

## write energy density figure 
pdf(paste0(write_figures, "Figure5.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
energy_density_plot
dev.off()


# plot momper 2017 MAG data -----------------------------------------------
source("momper2017.R")

## write MAG-level iron cycling bar plot figure 
pdf(paste0(write_figures, "Supp_Figure4.pdf"), height = 7.5, width = 13.33, pointsize = 1, paper = "a4r")
plot_grid(momper2017_taxa_function_plot, 
          plot_grid(momper2017_taxa_gene_feox, momper2017_taxa_gene_fered, align = "v", axis = "l", rel_widths = c(2, 2), ncol = 2), nrow = 2, rel_heights = c(2, 2.5))
dev.off()
