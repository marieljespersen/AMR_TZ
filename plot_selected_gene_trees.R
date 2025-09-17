library(treeio)
library(ggplot2)
library(ggtree)
library(stringr)
library(network)
library(RColorBrewer) 
library(gridExtra)

# set working dir
setwd("/Users/maloj/OneDrive - Danmarks Tekniske Universitet/HPC/ARG_flankophile/")

# metadata
flankophile_meta <- read.table(
  "1_hits_all.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)
#replace underscores in population with spaces
flankophile_meta$METADATA <- gsub("_", " ", flankophile_meta$METADATA)


### plot ####

# get gene list from directory
gene_files <- list.files("treefiles/")

custom_colors <- c(
  "Adult village population" = "#E74C3C",
  "Chimpanzee" = "#2E7F61",
  "Hadza population" = "#2E86C1",
  "School children population" = "#E67E22"
)

#f <- gene_files[1]
for(f in gene_files)
{
  g <- str_split(f, pattern = ".tree")[[1]][1]
  
  tree_file <- paste0("treefiles/", f)
  tree <- read.tree(tree_file)
  
  sub_meta <- flankophile_meta[which(flankophile_meta$OBSERVATION_ID %in% tree$tip.label),]
  Host_pop_tree <- rename_taxa(tree, sub_meta, OBSERVATION_ID, METADATA)
  Host <- c(Host_pop_tree$tip.label, rep(NA,Host_pop_tree$Nnode))
  x_lim <- 4*max(tree$edge.length)
  # regular tree
  p <- ggtree(tree)  + 
    geom_treescale(width=0.005) + 
    geom_tippoint(aes(color=Host), size=5) +
    ggtitle(g) + 
    theme(text = element_text(size = 20)) +
    theme(legend.position = "none") + 
    xlim(0, x_lim) +
    scale_color_manual(values = custom_colors)
  #print(p)
  plot_file <- paste0("results/tree_plots/tree_", g, ".plot.pdf")
  ggsave(plot=p,filename=plot_file,width=8, height=8, device="pdf", useDingbats=FALSE)
  print(g)
}

# end