#!/usr/bin/env Rscript

treefile = snakemake@input[[1]]

require(treeio) # Read in tree
require(ggtree) # Visualization
require(ggplot2)
require(ape) # Pairwise distances
require(pheatmap) # Distance heatmap
require(RColorBrewer)

tree <- read.newick(treefile)

tree_rect <- ggtree(tree, layout = 'rectangular', ladderize = F) +
  geom_treescale() +
  geom_tiplab(size=3, offset=0.0001, show.legend = F) +
  geom_nodelab(geom='label', size=2, fill='white', label.padding=unit(0.1, "lines"))

tree_circ <- ggtree(tree, layout = 'circular', ladderize = F) +
  geom_treescale() +
  geom_tiplab(size=3, offset=0.00025, show.legend = F)

tree_eqan <- ggtree(tree, layout = "equal_angle") +
  geom_tiplab(size=2, show.legend = F, hjust = -0.5) +
  coord_cartesian(clip="off")

tree_daylight <- ggtree(tree, layout = "daylight") +
  geom_tiplab(size=2, show.legend = F, hjust = -0.5) +
  coord_cartesian(clip="off")

outdir <- "figures/"
ggsave(paste0(outdir, "tree_rect.png"), tree_rect)
ggsave(paste0(outdir, "tree_circ.png"), tree_circ)
ggsave(paste0(outdir, "tree_eqan.png"), tree_eqan)
ggsave(paste0(outdir, "tree_daylight.png"), tree_daylight)

dists <- ape::cophenetic.phylo(as.phylo(tree))

pheatmap(dists, main = "Phylogenetic Tree Pairwise Distances",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         display_numbers = T, number_format = "%.4f",
         number_color = "black",
         clustering_method = "ward.D2",
         filename = paste0(outdir, "heatmap_distances.png"))
