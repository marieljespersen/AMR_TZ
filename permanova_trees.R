### analysis of trees compared to Host populations ####

# load libraries
library(tidyverse)
library(ape)
library(vegan)
library(ggExtra)
library(stringr)

# set working dir
setwd("/Users/maloj/OneDrive - Danmarks Tekniske Universitet/HPC/ARG_flankophile/")

# read metadata
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

### functions ####

#f <- files[[1]]
#  function to test tree
fit_adonis.f <- function(f) {
  
  print(f)
  # load tree
  
  tree<-read.tree(f)
 
  # get annotation
  md = flankophile_meta[match(tree$tip.label, flankophile_meta$OBSERVATION_ID),]
  colnames(md) = colnames(flankophile_meta)
  
  dup_Population <- unique(md$METADATA[duplicated(md$METADATA)])
  samples_to_drop <- md$OBSERVATION_ID[which(!md$METADATA %in% dup_Population)]
  
  # drop samples from where only one sample from the region is found
  tree <- drop.tip(tree, samples_to_drop)
  # create new sub metadata with only the samples left in the tree
  md <- md[!md$OBSERVATION_ID %in% samples_to_drop,]
  
  # check for multiple regions in tree before permanova
  if(length(unique(md$METADATA))>1)
  {
    # create  distance matrix from  tree
    PatristicDistMatrix<-cophenetic.phylo(tree)
    
    d = as.dist(PatristicDistMatrix)
    
    # testing distance according to conttinent (and longitude and latitude) using permanova
    fit = adonis2(d ~ METADATA, data=as.data.frame(md))
    return(fit)
  }
}


#### run tree files trough function ####
files <- list.files(path="treefiles", pattern = ".treefile", full.names=TRUE)
fit_list <- lapply(files, fit_adonis.f)
names(fit_list) <- lapply(files, function(f){strsplit(strsplit(f, '[/]')[[1]][2], '[.]')[[1]][1]})
fit_list <- compact(fit_list) # exclude NULL values from the trees that could not be tested (due tolacking multiple regions or one sample per region)


# extract p_vals and R2 values
pvals <- matrix(NA, nrow=length(fit_list), ncol=2)
for (i in 1:length(fit_list)) {
  pvals[i,] <- c(fit_list[[i]]$`Pr(>F)`[1], fit_list[[i]]$R2[1])
}
rownames(pvals) = names(fit_list)
colnames(pvals) <- c("p_val", "R2_model")

# adjust p values
pvals[,1] <- p.adjust(pvals[,1], method="BH")

pvals_df <- as.data.frame(pvals)
pvals_df$p_val <- round(pvals_df$p_val, digits = 3)
pvals_df$R2_model <- round(pvals_df$R2_model, digits = 3)

for(gene in rownames(pvals_df)){
  f <- paste0("treefiles/", gene, ".treefile")
  tree<-read.tree(f)
  n_samples <- length(tree$tip.label)
  
  #save to df
  pvals_df$n_samples[which(rownames(pvals_df)==gene)] <- n_samples
  
  print(gene)
  
}

write.table(pvals_df, file="results/permanova_res.txt", quote=FALSE)

# plot n_samples vs R2 coloured from significance level
pvals_df$significant <- ifelse(pvals_df$p_val < 0.05, "yes", "no")
pvals_df$gene <- rownames(pvals_df)
p <- ggplot(pvals_df, aes(x=n_samples, y=R2_model, color=significant, label=gene)) +
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values=c("yes"="darkgreen", "no"="steelblue")) +
  labs(title="PERMANOVA results: R2 vs number of samples per gene tree",
       x="Number of samples in gene tree",
       y="R2 of PERMANOVA model") +
  theme(text = element_text(size=15)) +
  geom_text(vjust=-1, size=3)
print(p)

# save plot
ggsave(filename="results/plots/permanova_R2_vs_n_samples.pdf", plot=p, width=8, height=6, device="pdf", useDingbats=FALSE)

