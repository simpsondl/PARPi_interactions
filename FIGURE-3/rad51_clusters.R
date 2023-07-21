library(readr)
library(WGCNA)
library(dplyr)
library(impute)
library(seriation)
library(ggraph)
library(tidygraph)
library(reshape2)
library(igraph)

# Load functions
source("helper_functions.R")

# Load data
gamma.gene.gis <- read_tsv(paste0("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/",
                                  "gene_combination_interaction_scores_gamma_oi_avg.txt"))
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
heatmap_clusters <- openxlsx::read.xlsx("../TABLES/Table_S4.xlsx")

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
                    minsupp = 4)

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


wg3 <- blockwiseModules(gene.mtx, power = 3, TOMType = "signed", minModuleSize = 3, numericLabels = TRUE)
wg7 <- blockwiseModules(gene.mtx, power = 7, TOMType = "signed", minModuleSize = 3, numericLabels = TRUE)

assigned.clusts <- data.frame(gene = rownames(gene.mtx), P3 = wg3$colors, P7 = wg7$colors)



rad51.fil <- heatmap_clusters$gene[heatmap_clusters$P6 == 1]
rad51.fil <- rad51.fil[!rad51.fil  %in% c("ZAR1L", "DCTN5")]
#filexp <- heatmap_clusters$gene[heatmap_clusters$Power4 == 2]
#filexp <- filexp[!filexp %in% c("ZAR1L", "DCTN5")]
#blm <- assigned.clusts$gene[assigned.clusts$P3 == 14]
blm <- heatmap_clusters$gene[heatmap_clusters$P4 == 10]
#spidr <- assigned.clusts$gene[assigned.clusts$P3 == 2]

big.rad51 <- assigned.clusts$gene[assigned.clusts$P3 == 2]
big.blm <- assigned.clusts$gene[assigned.clusts$P3 == 14]
big.rad51 <- big.rad51[!big.rad51  %in% c("ZAR1L", "DCTN5")]

# tmp.mtx.rad <- as.matrix(gene.mtx[rownames(gene.mtx) %in% rad51.fil, 
#                                           colnames(gene.mtx) %in% rad51.fil])
# # define an order for matrices
# roworder <- rownames(tmp.mtx.rad)
# colorder <- colnames(tmp.mtx.rad)
# cor.rows.mtx.rad <- 1-cor(t(tmp.mtx.rad))
# cor.cols.mtx.rad <- 1-cor(tmp.mtx.rad)
# d.rad <- as.dist(cor.rows.mtx.rad)
# dt.rad <- as.dist(cor.cols.mtx.rad)
# # Rearrange using OLO algorithm from seriation
# # These are the same when using OLO_average, but if you switch methods, they might not be
# o1.rad <- seriate(d.rad, method = "GW_average")
# o2.rad <- seriate(dt.rad, method = "GW_average")

# tmp.mtx.radexp <- as.matrix(gene.mtx[rownames(gene.mtx) %in% filexp, 
#                                   colnames(gene.mtx) %in% filexp])
# # define an order for matrices
# #roworder <- rownames(tmp.mtx.radexp)
# #colorder <- colnames(tmp.mtx.radexp)
# cor.rows.mtx.radexp <- 1-cor(t(tmp.mtx.radexp))
# cor.cols.mtx.radexp <- 1-cor(tmp.mtx.radexp)
# d.radexp <- as.dist(cor.rows.mtx.radexp)
# dt.radexp <- as.dist(cor.cols.mtx.radexp)
# # Rearrange using OLO algorithm from seriation
# # These are the same when using OLO_average, but if you switch methods, they might not be
# o1.radexp <- seriate(d.radexp, method = "OLO_average")
# o2.radexp <- seriate(dt.radexp, method = "OLO_average")


# tmp.mtx.blm <- as.matrix(gene.mtx[rownames(gene.mtx) %in% rad51.fil, 
#                                       colnames(gene.mtx) %in% blm])
# tmp.mtx.blm <- tmp.mtx.blm[roworder, ]
# cor.rows.mtx.blm <- 1-cor(t(tmp.mtx.blm))
# cor.cols.mtx.blm <- 1-cor(tmp.mtx.blm)
# d.blm <- as.dist(cor.rows.mtx.blm)
# dt.blm <- as.dist(cor.cols.mtx.blm)
# # Rearrange using OLO algorithm from seriation
# # These are the same when using OLO_average, but if you switch methods, they might not be
# o1.blm <- seriate(d.blm, method = "GW_average")
# o2.blm <- seriate(dt.blm, method = "GW_average")


tmp.mtx.rad.blm <- as.matrix(gene.mtx[rownames(gene.mtx) %in% c(rad51.fil, blm),
                                  colnames(gene.mtx) %in% c(rad51.fil, blm)])
# define an order for matrices
roworder.blm <- rownames(tmp.mtx.rad.blm)
colorder.blm <- colnames(tmp.mtx.rad.blm)
cor.rows.mtx.rad.blm <- 1-cor(t(tmp.mtx.rad.blm))
cor.cols.mtx.rad.blm <- 1-cor(tmp.mtx.rad.blm)
d.rad.blm <- as.dist(cor.rows.mtx.rad.blm)
dt.rad.blm <- as.dist(cor.cols.mtx.rad.blm)
# Rearrange using OLO algorithm from seriation
# These are the same when using OLO_average, but if you switch methods, they might not be
o1.rad.blm <- seriate(d.rad.blm, method = "GW_average")
o2.rad.blm <- seriate(dt.rad.blm, method = "GW_average")


tmp.mtx.big.rad.blm <- as.matrix(gene.mtx[rownames(gene.mtx) %in% c(big.rad51, big.blm),
                                      colnames(gene.mtx) %in% c(big.rad51, big.blm)])
# define an order for matrices
#roworder.blm <- rownames(tmp.mtx.rad.blm)
#colorder.blm <- colnames(tmp.mtx.rad.blm)
cor.rows.mtx.big.rad.blm <- 1-cor(t(tmp.mtx.big.rad.blm))
cor.cols.mtx.big.rad.blm <- 1-cor(tmp.mtx.big.rad.blm)
d.big.rad.blm <- as.dist(cor.rows.mtx.big.rad.blm)
dt.big.rad.blm <- as.dist(cor.cols.mtx.big.rad.blm)
# Rearrange using OLO algorithm from seriation
# These are the same when using OLO_average, but if you switch methods, they might not be
o1.big.rad.blm <- seriate(d.big.rad.blm, method = "GW_average")
o2.big.rad.blm <- seriate(dt.big.rad.blm, method = "GW_average")



tmp.mtx.rad <- as.matrix(gene.mtx[rownames(gene.mtx) %in% rad51.fil,
                                  colnames(gene.mtx) %in% rad51.fil])
#tmp.mtx.rad <- tmp.mtx.rad[roworder, ]
cor.rows.mtx.rad <- 1-cor(t(tmp.mtx.rad))
cor.cols.mtx.rad <- 1-cor(tmp.mtx.rad)
d.rad <- as.dist(cor.rows.mtx.rad)
dt.rad <- as.dist(cor.cols.mtx.rad)
# Rearrange using OLO algorithm from seriation
# These are the same when using OLO_average, but if you switch methods, they might not be
o1.rad <- seriate(d.rad, method = "GW_average")
o2.rad <- seriate(dt.rad, method = "OLO_average")

# Define color spectrum ranges
maxmag <- 2.723
col_fun <- colorRamp2(c(-1 * maxmag, 0, maxmag), c("dodgerblue", "white", "darkorange1"))

# Define annotation for the big map
tmp.anno <- assigned.clusts[assigned.clusts$gene %in% c(big.rad51, big.blm),1:2]
lvls <- c(rownames(tmp.mtx.big.rad.blm))
tmp.anno$gene <- factor(tmp.anno$gene, levels = lvls)
tmp.anno <- tmp.anno[order(tmp.anno$gene),]
tmp.anno$P3 <- ifelse(tmp.anno$P3 == 2, "RAD51 Filament", "Dissolvase Complex")

ha <- HeatmapAnnotation(Cluster = tmp.anno$P3, 
                        col = list(Cluster = c("RAD51 Filament" = "#004949", 
                                               "Dissolvase Complex" = "#6db6ff")),
                        show_annotation_name = FALSE)

# tmp.mtx.rad <- tmp.mtx.rad[o1.rad[[1]]$order,o1.rad[[1]]$order]
# tmp.mtx.blm <- tmp.mtx.blm[o1.rad[[1]]$order,o2.blm[[1]]$order]

# hrad <- Heatmap(tmp.mtx.rad,
#                   name = "GI",
#                   col = col_fun,
#                 cluster_rows = as.dendrogram(o1.rad[[1]]),
#                 cluster_columns = as.dendrogram(o1.rad[[1]]),
#                   show_row_dend = TRUE,
#                   show_column_dend = TRUE,
#                   #column_split = c(rep("A", 5), rep("B",12), rep("C",3)),
#                   row_names_side = "left",
#                   heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
#                   row_names_gp = gpar(fontsize = 7),
#                   column_names_gp = gpar(fontsize = 7),
#                   column_dend_height = unit(2, "cm"))

# hradexp <- Heatmap(tmp.mtx.radexp,
#                 name = "GI",
#                 col = col_fun,
#                 #cluster_rows = as.dendrogram(o1.radexp[[1]]), 
#                 #cluster_columns = as.dendrogram(o1.radexp[[1]]), 
#                 # row_order = 1:17,
#                 #column_order = 1:17,
#                 show_row_dend = TRUE,
#                 show_column_dend = TRUE,
#                 #column_split = c(rep("A", 5), rep("B",12), rep("C",3)),
#                 row_names_side = "left",
#                 heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
#                 row_names_gp = gpar(fontsize = 7),
#                 column_names_gp = gpar(fontsize = 7),
#                 column_dend_height = unit(2, "cm"))

# hblm <- Heatmap(tmp.mtx.blm,
#                 name = "GI",
#                 col = col_fun,
#                 row_order = 1:15,
#                 column_order = 1:5,
#                 show_row_dend = FALSE,
#                 show_column_dend = FALSE,
#                 #column_split = c(rep("A", 5), rep("B",12), rep("C",3)),
#                 heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
#                 row_names_side = "left",
#                 row_names_gp = gpar(fontsize = 7),
#                 column_names_gp = gpar(fontsize = 7),
#                 column_dend_height = unit(2, "cm"))
# 
# hlist <- hrad + hblm

# hrad <- Heatmap(tmp.mtx.rad,
#                 name = "GI",
#                 col = col_fun,
#                 cluster_rows = as.dendrogram(o1.rad[[1]]),
#                 cluster_columns = as.dendrogram(o1.rad[[1]]),
#                 #row_order = 1:17,
#                 #column_order = 1:6,
#                 show_row_dend = FALSE,
#                 show_column_dend = FALSE,
#                 #column_split = c(rep("A", 5), rep("B",12), rep("C",3)),
#                 heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
#                 row_names_side = "left",
#                 row_names_gp = gpar(fontsize = 7),
#                 column_names_gp = gpar(fontsize = 7),
#                 column_dend_height = unit(2, "cm"))

# def.ord <- hclust(d.rad, method = "average")
# tmp.mtx.blm <- tmp.mtx.blm[def.ord$order,]
# 
# hrad2 <- Heatmap(tmp.mtx.rad,
#                 name = "GI",
#                 col = col_fun,
#                 clustering_distance_columns = "pearson",
#                 clustering_distance_rows = "pearson",
#                 clustering_method_columns = "average",
#                 clustering_method_rows = "average",
#                 #cluster_rows = def.ord, 
#                 #cluster_columns = def.ord, 
#                 #row_order = 1:15,
#                 #column_order = 1:15,
#                 show_row_dend = FALSE,
#                 show_column_dend = TRUE,
#                 #column_split = c(rep("A", 5), rep("B",12), rep("C",3)),
#                 row_names_side = "left",
#                 heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
#                 row_names_gp = gpar(fontsize = 7),
#                 column_names_gp = gpar(fontsize = 7),
#                 column_dend_height = unit(2, "cm"))
# 
# 
# hblm2 <- Heatmap(tmp.mtx.blm,
#                 name = "GI",
#                 col = col_fun,
#                 clustering_distance_columns = "pearson",
#                 clustering_method_columns = "average",
#                 row_order = 1:15,
#                 #column_order = 1:5,
#                 show_row_dend = TRUE,
#                 show_column_dend = FALSE,
#                 #column_split = c(rep("A", 5), rep("B",12), rep("C",3)),
#                 heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
#                 row_dend_side = "right",
#                 row_names_side = "left",
#                 row_names_gp = gpar(fontsize = 7),
#                 column_names_gp = gpar(fontsize = 7),
#                 row_dend_width = unit(1.25, "cm"))
# 
# hlist2 <- hrad2 + hblm2


hradblm <- Heatmap(tmp.mtx.rad.blm,
                   rect_gp = gpar(type = "none"),
                name = "GI",
                col = col_fun,
                cluster_rows = as.dendrogram(o2.rad.blm[[1]]), 
                cluster_columns = as.dendrogram(o1.rad.blm[[1]]),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i == j){
                    grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                    # grid.text(rownames(tmp.mtx.rad.blm)[i], 
                    #           x = x, y = y,
                    #           gp = gpar(fontsize = 6))
                  } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                  } else {
                    if(abs(tmp.mtx.rad.blm[i,j]) >= 1.524){
                      #grid.rect(x,y,w,h,gp = gpar(col = "#bbbbbb"))
                      grid.circle(x = x, y = y,
                                  r = abs(tmp.mtx.rad.blm[i,j])/18 * min(unit.c(w, h)),
                                  gp = gpar(fill = col_fun(tmp.mtx.rad.blm[i,j]), col = NA))
                    }

                  }
                },
                show_row_dend = TRUE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                #row_dend_side = "right",
                heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                row_dend_side = "right",
                column_dend_height = unit(2, "cm"),
                show_heatmap_legend = FALSE)

hrad <- Heatmap(tmp.mtx.rad,
                   rect_gp = gpar(type = "none"),
                   name = "GI",
                   col = col_fun,
                   cluster_rows = as.dendrogram(o1.rad[[1]]), 
                   cluster_columns = as.dendrogram(o1.rad[[1]]),
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(i == j){
                       grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                     } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                     } else {
                       if(abs(tmp.mtx.rad[i,j]) >= 1.524){
                         grid.circle(x = x, y = y,
                                     r = abs(tmp.mtx.rad[i,j])/18 * min(unit.c(w, h)),
                                     gp = gpar(fill = col_fun(tmp.mtx.rad[i,j]), col = NA))
                       }
                       
                     }
                   },
                   show_row_dend = TRUE,
                   show_column_dend = FALSE,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   #row_dend_side = "right",
                   heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
                   row_names_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 7),
                   row_dend_side = "right",
                   column_dend_height = unit(2, "cm"),
                   show_heatmap_legend = FALSE)

hradlab <- Heatmap(tmp.mtx.rad,
                      rect_gp = gpar(type = "none"),
                      name = "GI",
                      col = col_fun,
                      cluster_rows = as.dendrogram(o2.rad[[1]]), 
                      cluster_columns = as.dendrogram(o2.rad[[1]]),
                      cell_fun = function(j, i, x, y, w, h, fill) {
                        if(i == j){
                          grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                          # grid.text(rownames(tmp.mtx.rad.blm)[i], 
                          #           x = x, y = y,
                          #           gp = gpar(fontsize = 6))
                        } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                          grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                        } else {
                          if(abs(tmp.mtx.rad[i,j]) >= 1.524){
                            #grid.rect(x,y,w,h,gp = gpar(col = "#bbbbbb"))
                            grid.circle(x = x, y = y,
                                        r = abs(tmp.mtx.rad[i,j])/18 * min(unit.c(w, h)),
                                        gp = gpar(fill = col_fun(tmp.mtx.rad[i,j]), col = NA))
                          }
                          
                        }
                      },
                      show_row_dend = TRUE,
                      show_column_dend = FALSE,
                      show_row_names = TRUE,
                      show_column_names = TRUE,
                      #row_dend_side = "right",
                      heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
                      row_names_gp = gpar(fontsize = 7),
                      column_names_gp = gpar(fontsize = 7),
                      column_names_side = "bottom",
                      row_dend_side = "right",
                      column_dend_height = unit(2, "cm"))

hradblmlab <- Heatmap(tmp.mtx.rad.blm,
                   rect_gp = gpar(type = "none"),
                   name = "GI",
                   col = col_fun,
                   cluster_rows = as.dendrogram(o2.rad.blm[[1]]), 
                   cluster_columns = as.dendrogram(o1.rad.blm[[1]]),
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(i == j){
                       grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                       # grid.text(rownames(tmp.mtx.rad.blm)[i], 
                       #           x = x, y = y,
                       #           gp = gpar(fontsize = 6))
                     } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                     } else {
                       if(abs(tmp.mtx.rad.blm[i,j]) >= 1.524){
                         #grid.rect(x,y,w,h,gp = gpar(col = "#bbbbbb"))
                         grid.circle(x = x, y = y,
                                     r = abs(tmp.mtx.rad.blm[i,j])/18 * min(unit.c(w, h)),
                                     gp = gpar(fill = col_fun(tmp.mtx.rad.blm[i,j]), col = NA))
                       }
                       
                     }
                   },
                   show_row_dend = TRUE,
                   show_column_dend = FALSE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   #row_dend_side = "right",
                   heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
                   row_names_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 7),
                   column_names_side = "bottom",
                   row_dend_side = "right",
                   column_dend_height = unit(2, "cm"))


hbigradblm <- Heatmap(tmp.mtx.big.rad.blm,
                   rect_gp = gpar(type = "none"),
                name = "GI",
                col = col_fun,
                cluster_rows = as.dendrogram(o2.big.rad.blm[[1]]), 
                cluster_columns = as.dendrogram(o1.big.rad.blm[[1]]),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i == j){
                    grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                    # grid.text(rownames(tmp.mtx.big.rad.blm)[i], 
                    #           x = x, y = y,
                    #           gp = gpar(fontsize = 6))
                  } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                  } else {
                    if(abs(tmp.mtx.big.rad.blm[i,j]) >= 1.5274){
                      grid.rect(x,y,w,h,gp = gpar(col = "#bbbbbb"))
                      grid.circle(x = x, y = y,
                                  r = abs(tmp.mtx.big.rad.blm[i,j])/19 * min(unit.c(w, h)),
                                  gp = gpar(fill = col_fun(tmp.mtx.big.rad.blm[i,j]), col = NA))
                    }

                  }
                },
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                top_annotation = ha,
                heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
                row_names_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 7),
                column_dend_height = unit(2, "cm"))

hbigradblmlab <- Heatmap(tmp.mtx.big.rad.blm,
                      rect_gp = gpar(type = "none"),
                      name = "GI",
                      col = col_fun,
                      cluster_rows = as.dendrogram(o2.big.rad.blm[[1]]), 
                      cluster_columns = as.dendrogram(o1.big.rad.blm[[1]]),
                      cell_fun = function(j, i, x, y, w, h, fill) {
                        if(i == j){
                          grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                          # grid.text(rownames(tmp.mtx.big.rad.blm)[i], 
                          #           x = x, y = y,
                          #           gp = gpar(fontsize = 6))
                        } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                          grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                        } else {
                          if(abs(tmp.mtx.big.rad.blm[i,j]) >= 1.524){
                            grid.rect(x,y,w,h,gp = gpar(col = "#bbbbbb"))
                            grid.circle(x = x, y = y,
                                        r = abs(tmp.mtx.big.rad.blm[i,j])/19 * min(unit.c(w, h)),
                                        gp = gpar(fill = col_fun(tmp.mtx.big.rad.blm[i,j]), col = NA))
                          }
                          
                        }
                      },
                      show_row_dend = FALSE,
                      show_column_dend = TRUE,
                      show_row_names = TRUE,
                      show_column_names = FALSE,
                      top_annotation = ha,
                      heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
                      row_names_gp = gpar(fontsize = 7),
                      column_names_gp = gpar(fontsize = 7),
                      column_dend_height = unit(2, "cm"))


# pdf("../FIGURES/FIGURE-3/radfil_blm.pdf", width = 8, height = 6)
# draw(hlist)
# dev.off()
# svg("../FIGURES/FIGURE-3/radfil_blm.svg", width = 8, height = 6)
# draw(hlist)
# dev.off()
# 
# pdf("../FIGURES/FIGURE-3/radfil_blm_noseriate.pdf", width = 8, height = 6)
# draw(hlist2)
# dev.off()
# svg("../FIGURES/FIGURE-3/radfil_blm_noseriate.svg", width = 8, height = 6)
# draw(hlist2)
# dev.off()

pdf("../FIGURES/FIGURE-3/radfil_blm_combinedclustsq.pdf", width = 6, height = 6)
draw(hradblm)
dev.off()
svg("../FIGURES/FIGURE-3/radfil_blm_combinedclustsq.svg", width = 6, height = 6)
draw(hradblm)
dev.off()

pdf("../FIGURES/FIGURE-3/radfil_blm_combinedclustsq_labels.pdf", width = 7, height = 6)
draw(hradblmlab)
dev.off()
svg("../FIGURES/FIGURE-3/radfil_blm_combinedclustsq_labels.svg", width = 7, height = 6)
draw(hradblmlab)
dev.off()

pdf("../FIGURES/FIGURE-S3/radfil_blm_expandedclustsq.pdf", width = 8, height = 7)
draw(hbigradblm)
dev.off()
svg("../FIGURES/FIGURE-S3/radfil_blm_expandededclustsq.svg", width = 8, height = 7)
draw(hbigradblm)
dev.off()

pdf("../FIGURES/FIGURE-S3/radfil_blm_expandedclustsq_labels.pdf", width = 8, height = 7)
draw(hbigradblmlab)
dev.off()
svg("../FIGURES/FIGURE-S3/radfil_blm_expandededclustsq_labels.svg", width = 8, height = 7)
draw(hbigradblmlab)
dev.off()

# pdf("../FIGURES/FIGURE-3/radexp.pdf", width = 7, height = 6)
# draw(hradexp)
# dev.off()
# svg("../FIGURES/FIGURE-3/radexp.svg", width = 7, height = 6)
# draw(hradexp)
# dev.off()
# 
pdf("../FIGURES/FIGURE-3/radfil_only.pdf", width = 6, height = 6)
draw(hrad)
dev.off()
svg("../FIGURES/FIGURE-3/radfil_blm.svg", width = 6, height = 6)
draw(hrad)
dev.off()
pdf("../FIGURES/FIGURE-3/radfil_only_labels.pdf", width = 5, height = 5)
draw(hradlab)
dev.off()
svg("../FIGURES/FIGURE-3/radfil_only_labels.svg", width = 5, height = 5)
draw(hradlab)
dev.off()



#################################################

mtx.melt <- melt(tmp.mtx.rad.blm)
mtx.melt$wt <- 0
mtx.melt$wt[mtx.melt$value > 1.524] <- 1
mtx.melt$wt[mtx.melt$value < -1.524] <- -1
mtx.melt$wt[mtx.melt$value < -2.723] <- -2
mtx.melt$wt[mtx.melt$value > 2.723] <- 2
mtx.melt$First <- apply(mtx.melt[,1:2], 1, min)
mtx.melt$Second <- apply(mtx.melt[,1:2], 1, max)
mtx.melt$withinclust <- !((mtx.melt$First %in% blm & mtx.melt$Second %in% rad51.fil) |
  (mtx.melt$First %in% rad51.fil & mtx.melt$Second %in% blm))
mtx.melt$inrad <- mtx.melt$First %in% rad51.fil & mtx.melt$Second %in% rad51.fil
mtx.melt$inblm <- mtx.melt$First %in% blm & mtx.melt$Second %in% blm
mtx.melt2 <- unique(mtx.melt[mtx.melt$wt != 0,c("First", "Second", "value", "wt", "withinclust", "inrad", "inblm")])

mtx.melt2$First[mtx.melt2$First == "NUDC" & mtx.melt2$Second == "PSMC3IP"] <- "PSMC3IP"
mtx.melt2$Second[mtx.melt2$First == "PSMC3IP" & mtx.melt2$Second == "PSMC3IP"] <- "NUDC"
mtx.melt2$First[mtx.melt2$First == "RAD51D" & mtx.melt2$Second == "SWI5"] <- "SWI5"
mtx.melt2$Second[mtx.melt2$First == "SWI5" & mtx.melt2$Second == "SWI5"] <- "RAD51D"
mtx.melt2$First[mtx.melt2$First == "BRCA2" & mtx.melt2$Second == "SWI5"] <- "SWI5"
mtx.melt2$Second[mtx.melt2$First == "SWI5" & mtx.melt2$Second == "SWI5"] <- "BRCA2"


gr <- as_tbl_graph(mtx.melt2) %>% 
  mutate(N = centrality_degree(mode = 'in')) 
gr <- gr %>%
  set_vertex_attr("type", value = vertex_attr(gr, "name") %in% blm)

# ggraph(gr, layout = "bipartite") + 
#   geom_edge_link(aes(width = as.character(abs(wt)), col = wt > 0, alpha = !withinclust)) + 
#   geom_edge_arc(aes(width = as.character(abs(wt)), col = wt > 0, alpha = inrad),
#                 fold = TRUE) + 
#   geom_node_point() + 
#   geom_node_text(aes(label = name), 
#                  repel = FALSE, 
#                  col = "black",
#                  size = 2) +
#   scale_edge_width_manual(values = c(.45, 1.25)) +
#   scale_edge_color_manual(values = c("dodgerblue", "darkorange1")) +
#   scale_color_manual(values = c("black", "#004949")) +
#   scale_edge_alpha_manual(values = c(0, 1)) +
#   theme_graph() +
#   theme(legend.position = "none")

g <- gr %>% 
  activate(nodes) %>%
  ggraph(layout = "graphopt") + 
  geom_edge_link(aes(width = as.character(abs(wt)), col = wt > 0)) + 
  geom_node_point(aes(col = !type), size = 9) + 
  geom_node_text(aes(label = name, col = type), 
                 repel = FALSE, 
                 size = 2, fontface = "bold") +
  scale_edge_width_manual(values = c(.45, 1.25)) +
  scale_edge_color_manual(values = c("dodgerblue", "darkorange1")) +
  scale_color_viridis(discrete = TRUE) +
  scale_edge_alpha_manual(values = c(0, 1)) +
  theme_graph() +
  theme(legend.position = "none")

l1 <- data.frame(x = g$data$x, y = g$data$y)

g2 <- gr %>%
  activate(nodes) %>%
  ggraph(layout = l1) + 
  geom_edge_link(aes(width = as.character(abs(wt)), col = wt > 0)) + 
  geom_node_point(aes(col = !type), size = 9) + 
  geom_node_text(aes(label = name, col = type), 
                 repel = FALSE, 
                 size = 2, fontface = "bold")

pdf("../FIGURES/FIGURE-3/rad51_blm_network.pdf", height = 5, width = 5)
g
dev.off()
