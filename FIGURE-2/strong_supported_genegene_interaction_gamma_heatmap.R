library(readr)
library(dplyr)
library(ggplot2)
library(impute)

# Load functions
source("helper_functions.R")

# Load data
gamma.gene.gis <- read_tsv(paste0("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/",
                                  "gene_combination_interaction_scores_gamma_oi_avg.txt"))
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")

# Get gene name from sgRNA.ID
single.pheno$Gene <- gsub("_.*", "", single.pheno$sgRNA.ID)
# Calculate single gene phenotypes
single.gene <- single.pheno %>% group_by(Gene) %>% summarise(Gamma = mean(Gamma.Avg))

# Merge in gene names
gamma.gene.gis <- inner_join(gamma.gene.gis,
                             id.map)

# Identify strong, supported interactions
ssi <- select_genes(scores = gamma.gene.gis,
                    minmagnitude = 2.723,
                    minstrong = 4,
                    minsupp = 3)

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol.gamma.gene.gis <- gamma.gene.gis[!grepl("non-targeting", gamma.gene.gis$GeneCombinationName),]

# Make a 538 x 538 grid with gene names
gene.grid <- expand.grid(unique(noncontrol.gamma.gene.gis$Gene1), unique(noncontrol.gamma.gene.gis$Gene1))
# Merge in interaction scores
gene.grid <- left_join(gene.grid, 
                        noncontrol.gamma.gene.gis, 
                        by = c("Var1" = "Gene2", "Var2" = "Gene1"))
gene.grid <- left_join(gene.grid, 
                       noncontrol.gamma.gene.gis, 
                       by = c("Var1" = "Gene1", "Var2" = "Gene2"))
# Combine the columns
gene.grid$Gene.GI <- ifelse(is.na(gene.grid$InteractionScore.x), gene.grid$InteractionScore.y, gene.grid$InteractionScore.x)
gene.grid$N <- ifelse(is.na(gene.grid$N.x), gene.grid$N.y, gene.grid$N.x)

# Remove extraneous columns
gene.grid <- gene.grid[,c("Var1", "Var2", "Gene.GI", "N")]
gene.grid$Gene.GI <- as.numeric(gene.grid$Gene.GI)
gene.grid$N <- as.numeric(gene.grid$N)

# Go from long to wide
gene.mtx <- pivot_wider(gene.grid[,1:3], names_from = Var2, values_from = Gene.GI)
gene.mtx <- as.data.frame(gene.mtx)
### Rename rows
rownames(gene.mtx) <- gene.mtx$Var1
### Remove redundant column
gene.mtx <- gene.mtx[,2:ncol(gene.mtx)]

# Impute missing values
gene.mtx <- impute.knn(as.matrix(gene.mtx))$data

# Make triangular heatmaps
gamma.ssi.tri <- create_lower_triangle_hm(mtx = gene.mtx,
                                    maxmag = 3.5,
                                    genes = ssi,
                                    annos = single.gene,
                                    include_anno = TRUE)

# Save heatmaps
pdf("../FIGURES/FIGURE-2/figure_2b_gamma_ssi_triangular_no_labels.pdf", width = 10, height = 10)
draw(gamma.ssi.tri[[2]])
dev.off()
svg("../FIGURES/FIGURE-2/figure_2b_gamma_ssi_triangular_no_labels.svg", width = 10, height = 10)
draw(gamma.ssi.tri[[2]])
dev.off()
pdf("../FIGURES/FIGURE-2/figure_2b_gamma_ssi_triangular_labels.pdf", width = 10, height = 10)
draw(gamma.ssi.tri[[1]])
dev.off()
svg("../FIGURES/FIGURE-2/figure_2b_gamma_ssi_triangular_labels.svg", width = 10, height = 10)
draw(gamma.ssi.tri[[1]])
dev.off()

