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

tmp.mtx <- gene.mtx[ssi, ssi]
d <- as.dist(1 - cor(as.matrix(tmp.mtx)))
dt <- as.dist(t(1 - cor(as.matrix(tmp.mtx))))
# Rearrange using OLO algorithm from seriation
o1 <- seriate(d, method = "OLO_average")
o2 <- seriate(dt, method = "GW_average")
maxmag <- 3.5
# Define color spectrum ranges
col_fun <- colorRamp2(c(-1 * maxmag, 0, maxmag), c("dodgerblue", "white", "darkorange1"))

sg4 <- blockwiseModules(tmp.mtx, power = 4, TOMType = "signed", minModuleSize = 3, numericLabels = TRUE)
sg6 <- blockwiseModules(tmp.mtx, power = 6, TOMType = "signed", minModuleSize = 3, numericLabels = TRUE)

clusters <- data.frame(gene = rownames(tmp.mtx), P4 = sg4$colors, P6 = sg6$colors)

clusters$C4 <- clusters$P4
clusters$C4[clusters$C4 == 0] <- NA
clusters$C6 <- clusters$P6
clusters$C6[clusters$C6 == 0] <- NA
clusters$gene <- factor(clusters$gene, levels = ssi)
clusters <- clusters[order(clusters$gene),]
clusters$C4 <- factor(as.character(clusters$C4), levels = as.character(1:9))
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
                                   C4 = setNames(paletteMartin[c(2,4,7,3,11,14,13,10,5,8,6,9,12,15)], as.character(1:14))), 
                        na_col = "white")
#col
c.both2 <- HeatmapAnnotation(Gamma = tmp.anno$Gamma,
                             C4 = clusters$C4,
                             C6 = clusters$C6,
                             show_annotation_name = FALSE,
                             show_legend = FALSE,
                             col = list(Gamma = colorRamp2(c(-.45,.1), c("darkred", "white")),
                                        C6 = setNames(paletteMartin[c(14,13,6,7,11,4,2)], as.character(7:1)), 
                                        C4 = setNames(paletteMartin[c(2,4,7,3,11,14,13,10,5,8,6,9,12,15)], as.character(1:14))), 
                             na_col = "white")

c3 <- rowAnnotation(C4 = clusters$C4,
                    show_annotation_name = FALSE, 
                    show_legend = FALSE,
                    col = list(C4 = setNames(paletteMartin[c(2,4,7,3,11,14,13,10,5,8,6,9,12,15)], as.character(1:14))), 
                    na_col = "white")

hm <- Heatmap(tmp.mtx, 
        rect_gp = gpar(type = "none"), 
        cluster_rows = as.dendrogram(o2[[1]]), 
        cluster_columns = as.dendrogram(o2[[1]]), 
        show_row_dend = TRUE,
        show_column_dend = FALSE,
        row_dend_side = "left",
        heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i == j){
            grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
          } else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        },
        col = col_fun,
        name = "GI",
        bottom_annotation = c.both2, 
        left_annotation = c3,
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 4),
        column_dend_height = unit(2, "cm"))

hm.nolab <- Heatmap(tmp.mtx, 
              rect_gp = gpar(type = "none"), 
              cluster_rows = as.dendrogram(o2[[1]]), 
              cluster_columns = as.dendrogram(o2[[1]]), 
              show_row_dend = TRUE,
              show_column_dend = FALSE,
              row_dend_side = "left",
              heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i == j){
                  grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                } else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
              },
              col = col_fun,
              name = "GI",
              bottom_annotation = c.both2, 
              left_annotation = c3,
              column_names_side = "bottom",
              show_row_names = FALSE,
              show_column_names = FALSE,
              show_heatmap_legend = FALSE,
              row_dend_width = unit(2, "cm"))

pdf("../FIGURES/FIGURE-2/heatmap_with_cluster_annotations.pdf", width = 10, height = 10)
draw(hm)
dev.off()

svg("../FIGURES/FIGURE-2/heatmap_with_cluster_annotations.svg", width = 10, height = 10)
draw(hm)
dev.off()

pdf("../FIGURES/FIGURE-2/heatmap_with_cluster_annotations_nolabels.pdf", width = 10, height = 10)
draw(hm.nolab)
dev.off()

svg("../FIGURES/FIGURE-2/heatmap_with_cluster_annotations_nolabels.svg", width = 10, height = 10)
draw(hm.nolab)
dev.off()
