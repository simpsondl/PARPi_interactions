library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(impute)
library(ComplexHeatmap)
library(circlize)
library(seriation)

#tauclusts <- read_tsv("../DATA/tau_clusters_p4p13.txt")

# Defined by clustering all genes on tau interaction scores
#fancs <- tauclusts$gene[tauclusts$P13 == 2]
#brca1a <- c("BABAM1", "BRCA1", "BRCC3", "FAM175A", "UIMC1")
#buffs <- tauclusts$gene[tauclusts$P4 == 7]

# DEFINED BY LITERATURE REVIEW
fancs <- c("FANCB", "FANCA", "FANCG", "FANCF", "FANCM", "C19orf40",  
           "FANCE", "FANCC", "FANCD2", "FANCI", "UBE2T",
           "BRIP1", "RAD18")
brca1a <- c("UIMC1", "BRCC3", "BABAM1", "FAM175A", "BRCA1")
buffs <- c("TP53BP1", "RIF1", "MAD2L2", "FAM35A", "C20orf196", "CCAR1",  "CTC1", 
            "OBFC1", "TEN1", "ATM")

#ssi <- c(unique(c(fancs, brca1a, buffs)), "AUNIP")

tau.gene.gis <- read_tsv(paste0("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/",
                                  "gene_combination_interaction_scores_tau_oi_avg.txt"))
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
rho.single <- read_tsv("../DATA/single_sgRNA_rho_phenotypes.txt")

tau.gene.gis <- inner_join(tau.gene.gis,
                             id.map)
rho.single$gene <- gsub("_.*", "", rho.single$sgRNA.ID)

rho.gene <- rho.single %>% group_by(gene) %>%
  summarise(rho = mean(Rho.Avg))

noncontrol.tau.gene.gis <- tau.gene.gis[!grepl("non-targeting", tau.gene.gis$GeneCombinationName),]

# Make a 538 x 538 grid with gene names
gene.grid <- expand.grid(unique(noncontrol.tau.gene.gis$Gene1), unique(noncontrol.tau.gene.gis$Gene1))
# Merge in interaction scores
gene.grid <- left_join(gene.grid, 
                       noncontrol.tau.gene.gis, 
                       by = c("Var1" = "Gene2", "Var2" = "Gene1"))
gene.grid <- left_join(gene.grid, 
                       noncontrol.tau.gene.gis, 
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

tmp.mtx <- gene.mtx[c("AUNIP", "BRCA1", fancs), 
                    c("AUNIP", c(buffs, brca1a, fancs))]

d <- as.dist(1 - cor(as.matrix(t(tmp.mtx))))
dt <- as.dist(1 - cor(as.matrix(tmp.mtx)))
# Rearrange using OLO algorithm from seriation
o1 <- seriate(d, method = "OLO_average")
o2 <- seriate(dt, method = "OLO_average")
maxmag <- 3.5
# Define color spectrum ranges
col_fun <- colorRamp2(c(-1 * maxmag, 0, maxmag), c("dodgerblue", "white", "darkorange1"))

rho.sub <- rho.gene[rho.gene$gene %in% c("AUNIP", c(buffs, brca1a)),]
rho.sub$gene <- factor(rho.sub$gene, levels = c("AUNIP", c(buffs, brca1a)))
rho.sub <- rho.sub[order(rho.sub$gene),]
rho.sub$col <- "black"
rho.sub$col[rho.sub$gene %in% buffs] <- "#CC6677" 
rho.sub$col[rho.sub$gene %in% brca1a] <- "#DDCC77"
#ha <- HeatmapAnnotation(Rho = anno_barplot(rho.sub$rho,
#                                           bar_width = .9,
#                                           gp = gpar(fill = rho.sub$col,
#                                                     col = rho.sub$col)))

tmp.mtx <- tmp.mtx[,c("AUNIP", 
                      "BRCA1", "FAM175A", "UIMC1", "BRCC3", "BABAM1",
                      "RIF1", "ATM", "MAD2L2", "TEN1",
                      "OBFC1", "FAM35A", "C20orf196", "CTC1", 
                      "CCAR1", "TP53BP1",
                      'FANCB', 'FANCA', 'C19orf40', 'FANCD2', 'FANCI', 'UBE2T', 
                      'FANCM', 'FANCF', 'FANCG', 'FANCE', 'RAD18', 'FANCC', 'BRIP1')]

tmp.mtx2 <- tmp.mtx
tmp.mtx2["AUNIP", "AUNIP"] <- NA
tmp.mtx2["BRCA1", "BRCA1"] <- NA
tmp.mtx2["FANCA", "FANCA"] <- NA
tmp.mtx2["FANCB", "FANCB"] <- NA
tmp.mtx2["FANCC", "FANCC"] <- NA
tmp.mtx2["FANCD2", "FANCD2"] <- NA
tmp.mtx2["FANCE", "FANCE"] <- NA
tmp.mtx2["FANCF", "FANCF"] <- NA
tmp.mtx2["FANCG", "FANCG"] <- NA
tmp.mtx2["FANCI", "FANCI"] <- NA
tmp.mtx2["BRIP1", "BRIP1"] <- NA
tmp.mtx2["RAD18", "RAD18"] <- NA
tmp.mtx2["FANCM", "FANCM"] <- NA
tmp.mtx2["C19orf40", "C19orf40"] <- NA
tmp.mtx2["UBE2T", "UBE2T"] <- NA

hm <- Heatmap(tmp.mtx2,  
        #cluster_rows = as.dendrogram(o1[[1]]), 
        #cluster_columns = as.dendrogram(o2[[1]]),
        column_order = c("AUNIP", 
                         "BRCA1", "FAM175A", "UIMC1", "BRCC3", "BABAM1",
                         "RIF1", "ATM", "MAD2L2", "TEN1",
                         "OBFC1", "FAM35A", "C20orf196", "CTC1", 
                         "CCAR1", "TP53BP1",
                         'FANCB', 'FANCA', 'C19orf40', 'FANCD2', 'FANCI', 'UBE2T', 
                         'FANCM', 'FANCF', 'FANCG', 'FANCE', 'RAD18', 'FANCC', 'BRIP1'),
        row_order = c("AUNIP",
                      "BRCA1",
                      'FANCB', 'FANCA', 'C19orf40', 'FANCD2', 'FANCI', 'UBE2T', 
                      'FANCM', 'FANCF', 'FANCG', 'FANCE', 'RAD18', 'FANCC', 'BRIP1'
                      ),
        row_dend_side = "right", 
        #column_dend_side = "top",
        heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
        row_names_side = "right",
        column_split = c("A", rep("B", 5), rep("C", 10), rep ("D", 13)),
        row_split = c("A", "B", rep("C", 13)),
        column_gap = unit(3, "mm"),
        col = col_fun,
        name = "GI",
        column_names_side = "bottom",
        #bottom_annotation = ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))

pdf("../FIGURES/FIGURE-6/aunip_v_resection_hm_072023.pdf", width = 6, height = 3)
draw(hm)
dev.off()

##########################################

# Compare clustering if distances from full dataset are used
dfull <- as.dist(1 - cor(as.matrix(gene.mtx)))
dsub1 <- as.matrix(dfull)[c("AUNIP",c(buffs, brca1a)), c("AUNIP",c(buffs, brca1a))]
dsub2 <- as.matrix(dfull)[c("AUNIP", fancs), c("AUNIP", fancs)]
dsub1 <- as.dist(dsub1)
dsub2 <- as.dist(dsub2)

osub1 <- seriate(dsub1, method = "OLO_average")
osub2 <- seriate(dsub2, method = "OLO_average")
