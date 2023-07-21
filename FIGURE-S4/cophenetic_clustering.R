library(readr)
library(dplyr)
library(ggplot2)
library(impute)
library(seriation)
library(dendextend)

# Load functions
source("helper_functions.R")

# Load data
gamma.gene.gis <- read_tsv(paste0("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/",
                                  "gene_combination_interaction_scores_gamma_oi_avg.txt"))
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
clusts <- read_tsv("../DATA/clusts_0717.txt")

# Merge in gene names
gamma.gene.gis <- inner_join(gamma.gene.gis,
                             id.map)

# Identify strong, supported interactions
# ssi <- select_genes(scores = gamma.gene.gis,
#                     minmagnitude = 2.723,
#                     minstrong = 4,
#                     minsupp = 4)

tmp <- gamma.gene.gis[(gamma.gene.gis$InteractionScore >= 2.7555 | gamma.gene.gis$InteractionScore <= -2.754) &
                        gamma.gene.gis$N >= 4,] %>%
  filter(Category == "X+Y")
gene.list <- unique(c(tmp$Gene1, tmp$Gene2))
sg.df <- data.frame(gene = gene.list, freq = NA)
for(i in sg.df$gene){
  sg.df$freq[sg.df$gene == i] <- sum(tmp$Gene1 == i | tmp$Gene2 == i)
}

ssi <- sg.df$gene[sg.df$freq >= 4]

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

d <- as.dist(1 - cor(as.matrix(gene.mtx)))
dsub <- as.matrix(d)[ssi,ssi]
dsub <- as.dist(dsub)

tmp.mtx <- gene.mtx[ssi, ssi]
dtmp <- as.dist(1 - cor(as.matrix(tmp.mtx)))

o1 <- seriate(dtmp, method = "OLO_average")
o2 <- seriate(dsub, method = "OLO_average")

clusts$C4 <- clusts$P4
clusts$C4[clusts$C4 == 0] <- NA
clusts$gene <- factor(clusts$gene, levels = ssi[o1[[1]]$order])
clusts <- clusts[order(clusts$gene),]
clusts$C4 <- factor(as.character(clusts$C4), levels = as.character(1:13))
clusts$col <- sapply(clusts$P4, function(x) c("white",paletteMartin[c(2,4,7,3,11,14,13,10,5,8,6,9,12,15)])[x+1])

clusts2 <- clusts
clusts2$gene <- factor(clusts2$gene, levels = ssi[o2[[1]]$order])
clusts2 <- clusts2[order(clusts2$gene),]

clusts3 <- clusts
clusts3$gene <- factor(clusts3$gene, levels = ssi)
clusts3 <- clusts3[order(clusts3$gene),]

# an <- rowAnnotation(SUB = clusts$C4,
#                     FULL = clusts2$C4,
#                     col = list(SUB = setNames(paletteMartin[c(3,2,4,7,11,6,5,13,8,9,10,14,12,15)], 
#                                              as.character(1:14)),
#                                FULL = setNames(paletteMartin[c(3,2,4,7,11,6,5,13,8,9,10,14,12,15)], 
#                                                as.character(1:14))), 
#                     na_col = "#bbbbbb")
# draw(an)

temp_col <- clusts$col
temp_col <- factor(temp_col, unique(temp_col))

den <- as.dendrogram(o1[[1]])
den <- color_branches(den, clusters = clusts$P4, 
                      col = c(paletteMartin[c(4,2,5,10,7,14,11,3,13,8,6,9,12,15)],"white"))

svg("../FIGURES/FIGURE-S4/cophenetic_full_vs_sub.svg", height = 4, width = 6)
den %>%
  set_labels(NULL) %>%
  raise.dendrogram(heiget_to_add = -.15) %>%
  plot(axes = FALSE)
colored_bars(colors = cbind(clusts2$col, clusts$col), 
             dend = den, 
             sort_by_labels_order = FALSE,
             rowLabels = c("FULL", "SUB"), 
             y_shift = -.02)
dev.off()

cor_cophenetic(as.dendrogram(o1[[1]]), as.dendrogram(o2[[1]]))
