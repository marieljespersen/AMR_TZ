
library(stringr)
library(ggplot2)
library(dplyr)

# set working directory
setwd("/Users/maloj/OneDrive - Danmarks Tekniske Universitet/HPC/ARG_flankophile/")

#### read data ####
flanko_res <- read.table("1_hits_all.tsv", header=TRUE)
colnames(flanko_res) <- c("FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS", "COVERAGE_P", "IDENTITY_P", "CONTIG_LENGTH", "Population", "ASSEMBLY_NAME", "OBSERVATION_ID", "VARIANT")

# Clean up
flanko_res$Population <- str_replace_all(flanko_res$Population, "_", " ")
flanko_res <- flanko_res[!(flanko_res$Population %in% c("Antelopes", "mislabeled")), ]
flanko_res$GENE <- str_replace_all(flanko_res$GENE, "gene_", "")

# Filter genes with at least 10 observations
gene_counts <- flanko_res %>%
  group_by(GENE) %>%
  summarise(n = n()) %>%
  filter(n >= 10)

filtered_genes <- gene_counts$GENE
flanko_res <- flanko_res[flanko_res$GENE %in% filtered_genes, ]

# Custom colors
custom_colors <- c(
  "Adult village population" = "#E74C3C",
  "Chimpanzee" = "#2E7F61",
  "Hadza population" = "#2E86C1",
  "School children population" = "#E67E22"
)

# Create PDF with one plot per page
pdf("results/plots/gene_identity_per_gene_min10.pdf", width=6, height=4)

for (gene in unique(flanko_res$GENE)) {
  gene_data <- flanko_res[flanko_res$GENE == gene, ]
  total_n <- nrow(gene_data)
  
  # Count per population
  pop_counts <- gene_data %>%
    group_by(Population) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(Population, "\n(n=", n, ")"))
  
  # Merge labels back into data
  gene_data <- merge(gene_data, pop_counts[, c("Population", "label")], by = "Population")
  
  p <- ggplot(gene_data, aes(x = label, y = IDENTITY_P, fill = Population)) +
    geom_boxplot() +
    labs(title = paste("Gene:", gene, "(n =", total_n, ")")) +
    theme_bw() +
    ylab("% Identity") +
    xlab("Population") +
    ylim(99, 100) +
    scale_fill_manual(values = custom_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          text = element_text(size = 10))
  
  print(p)
}

dev.off()
