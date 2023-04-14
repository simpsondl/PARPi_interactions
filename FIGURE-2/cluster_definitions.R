library(readr)
library(WGCNA)
library(dplyr)
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

# Cluster data and cut tree dynamically
tmp.mtx <- as.matrix(gene.mtx[rownames(gene.mtx) %in% ssi, colnames(gene.mtx) %in% ssi])
d <- as.dist(1 - cor(as.matrix(tmp.mtx)))
o1 <- seriate(d, method = "OLO_average")

# geneclust <- hclust(d, method = "average")
# dyn <- cutreeDynamic(dendro = geneclust, 
#                      distM = as.matrix(d), 
#                      deepSplit = 2, 
#                      pamRespectsDendro = TRUE, 
#                      minClusterSize = 3,
#                      method = "hybrid", 
#                      pamStage = TRUE)

#WGCNA
wg3 <- blockwiseModules(tmp.mtx, 
                        power = 3, 
                        TOMType = "unsigned",
                        minModuleSize = 3, 
                        numericLabels = TRUE) 

wg4 <- blockwiseModules(tmp.mtx, 
                        power = 4, 
                        TOMType = "signed",
                        minModuleSize = 3,
                        #deepSplit = 4,
                        numericLabels = TRUE) 

wg6 <- blockwiseModules(tmp.mtx, 
                       power = 6, 
                       TOMType = "unsigned",
                       minModuleSize = 3, 
                       numericLabels = TRUE) 

assigned.annotations <- data.frame(gene = rownames(tmp.mtx),
                                   P3 = wg3$colors,
                                   P4 = wg4$colors,
                                   P6 = wg6$colors)

for_labeling <- data.frame(gene = assigned.annotations$gene, cluster = NA)
for_labeling$cluster[assigned.annotations$P6 != 0] <- assigned.annotations$P6[assigned.annotations$P6 != 0]
for_labeling$cluster[assigned.annotations$P4 == 0] <- 0
for_labeling$cluster[assigned.annotations$P4 == 14] <- 14
for_labeling$cluster[assigned.annotations$P4 == 13] <- 13
for_labeling$cluster[assigned.annotations$P4 == 11] <- 12
for_labeling$cluster[assigned.annotations$P4 == 10] <- 11
for_labeling$cluster[assigned.annotations$P4 == 9] <- 10
for_labeling$cluster[assigned.annotations$P4 == 7] <- 9
for_labeling$cluster[assigned.annotations$P4 == 1] <- 8
for_labeling$cluster[is.na(for_labeling$cluster)] <- 0

write_tsv(for_labeling, "../DATA/heatmap_clusters_power6base.txt")

for_labeling <- data.frame(gene = assigned.annotations$gene, cluster = NA)
for_labeling$cluster[assigned.annotations$P6 != 0] <- assigned.annotations$P6[assigned.annotations$P6 != 0]
for_labeling$cluster[assigned.annotations$P4 == 0] <- 0
for_labeling$cluster[assigned.annotations$P4 == 14] <- 14
for_labeling$cluster[assigned.annotations$P4 == 13] <- 13
for_labeling$cluster[assigned.annotations$P4 == 11] <- 12
for_labeling$cluster[assigned.annotations$P4 == 10] <- 11
for_labeling$cluster[assigned.annotations$P4 == 9] <- 10
for_labeling$cluster[assigned.annotations$P4 == 7] <- 9
for_labeling$cluster[assigned.annotations$P4 == 1] <- 8

for_labeling$cluster[assigned.annotations$P4 == 12] <- 7
for_labeling$cluster[assigned.annotations$P4 == 8] <- 6
for_labeling$cluster[assigned.annotations$P4 == 6] <- 5
for_labeling$cluster[is.na(for_labeling$cluster)] <- 0

write_tsv(for_labeling, "../DATA/heatmap_clusters_additionalrescue.txt")

