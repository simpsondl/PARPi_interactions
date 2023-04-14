library(readr)
library(dplyr)
library(ggplot2)
library(impute)
library(colorBlindness)

# Load functions
source("helper_functions.R")

# Load data
gamma.gene.gis <- read_tsv(paste0("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/",
                                  "gene_combination_interaction_scores_gamma_oi_avg.txt"))
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")

clusters <- read_delim("C:/Users/danny/Desktop/GI_PUBLICATION/DATA/heatmap_clusters_power4power6.txt",
delim = "\t", escape_double = FALSE,
trim_ws = TRUE)


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

tmp.mtx <- gene.mtx[ssi, ssi]
d <- as.dist(1 - cor(as.matrix(tmp.mtx)))
dt <- as.dist(t(1 - cor(as.matrix(tmp.mtx))))
# Rearrange using OLO algorithm from seriation
o1 <- seriate(d, method = "OLO_average")
o2 <- seriate(dt, method = "GW_average")
maxmag <- 3.5
# Define color spectrum ranges
col_fun <- colorRamp2(c(-1 * maxmag, 0, maxmag), c("dodgerblue", "white", "darkorange1"))

clusters$C4 <- clusters$Power4
clusters$C4[clusters$C4 == 0] <- NA
clusters$C6 <- clusters$Power6
clusters$C6[clusters$C6 == 0] <- NA
clusters$gene <- factor(clusters$gene, levels = ssi)
clusters <- clusters[order(clusters$gene),]
clusters$C4 <- factor(as.character(clusters$C4), levels = as.character(1:14))
clusters$C6 <- factor(as.character(clusters$C6), levels = as.character(1:7))


tmp.anno <- single.gene[single.gene$Gene %in% ssi,]
tmp.anno$Gene <- factor(tmp.anno$Gene, levels = ssi)
tmp.anno <- tmp.anno[order(tmp.anno$Gene),]

#row
c.both <- rowAnnotation(C6 = clusters$C6,
                        C4 = clusters$C4,
                        Gamma = tmp.anno$Gamma,
                        show_annotation_name = FALSE, 
                        col = list(Gamma = colorRamp2(c(-.45,.1), c("darkred", "white")),
                                   C6 = setNames(paletteMartin[c(14,13,6,7,11,4,2)], as.character(7:1)), 
                                   C4 = setNames(paletteMartin[c(3,2,4,7,11,6,5,13,8,9,10,14,12,15)], as.character(1:14))), 
                        na_col = "white")
#col
c.both2 <- HeatmapAnnotation(C4 = as.character(clusters$C4), 
                             C6 = as.character(clusters$C6), 
                             show_annotation_name = FALSE, 
                             #col = list(C6 = test.names, C4 = test.names2), 
                             na_col = "white")

hm <- Heatmap(tmp.mtx, 
        rect_gp = gpar(type = "none"), 
        cluster_rows = as.dendrogram(o1[[1]]), 
        cluster_columns = as.dendrogram(o1[[1]]), 
        show_row_dend = FALSE, 
        column_dend_side = "bottom",
        heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        },
        col = col_fun,
        name = "GI",
        left_annotation = c.both, 
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 4),
        column_dend_height = unit(2, "cm"))

hm.nolab <- Heatmap(tmp.mtx, 
              rect_gp = gpar(type = "none"), 
              cluster_rows = as.dendrogram(o1[[1]]), 
              cluster_columns = as.dendrogram(o1[[1]]), 
              show_row_dend = FALSE, 
              column_dend_side = "bottom",
              heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
              },
              col = col_fun,
              name = "GI",
              left_annotation = c.both, 
              column_names_side = "bottom",
              show_row_names = FALSE,
              show_column_names = FALSE,
              column_dend_height = unit(2, "cm"))

pdf("../FIGURES/FIGURE-2/figure_2b_cluster_annotations.pdf", width = 10, height = 10)
draw(hm)
dev.off()

svg("../FIGURES/FIGURE-2/figure_2b_cluster_annotations.svg", width = 10, height = 10)
draw(hm)
dev.off()

pdf("../FIGURES/FIGURE-2/figure_2b_cluster_annotations_nolabels.pdf", width = 10, height = 10)
draw(hm.nolab)
dev.off()

svg("../FIGURES/FIGURE-2/figure_2b_cluster_annotations_nolabels.svg", width = 10, height = 10)
draw(hm.nolab)
dev.off()
