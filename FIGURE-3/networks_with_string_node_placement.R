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
string_placement <- read_delim("../DATA/External_Data/string_network_coordinates (4).tsv",
                               delim = "\t", escape_double = FALSE, trim_ws = TRUE)
string_links <- read_table("../DATA/External_Data/9606.protein.links.full.v11.5.txt")
string_interactions_short <- read.delim("../DATA/External_Data/string_interactions_short (4).tsv")
heatmap_clusters <- read_tsv("../DATA/clusts_0717.txt")


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


rad51.fil <- heatmap_clusters$gene[heatmap_clusters$P6 == 1]
rad51.fil <- rad51.fil[!rad51.fil  %in% c("ZAR1L", "DCTN5")]
#filexp <- heatmap_clusters$gene[heatmap_clusters$Power4 == 2]
#filexp <- filexp[!filexp %in% c("ZAR1L", "DCTN5")]
#blm <- assigned.clusts$gene[assigned.clusts$P3 == 14]
blm <- heatmap_clusters$gene[heatmap_clusters$P4 == 8]

# Get node placement from STRING
string_placement <- string_placement[,c(1,3,4)]
colnames(string_placement) <- c("gene", "x", "y")
string_placement$gene[string_placement$gene == "NSMCE3"] <- "NDNL2"
string_placement <- string_placement[string_placement$gene %in% c(rad51.fil, blm),]
string_placement$y[string_placement$gene == "TOP3A"] <- .63
string_placement$y[string_placement$gene == "NDNL2"] <- .54
string_placement$x[string_placement$gene == "XRCC2"] <- .45
string_placement$y[string_placement$gene == "XRCC2"] <- .7
string_placement$y[string_placement$gene == "RAD51B"] <- .5
string_placement$x[string_placement$gene == "NUDC"] <- .42
string_placement$x[string_placement$gene == "RAD51D"] <- .49
string_placement$y[string_placement$gene == "RAD51D"] <- .605
string_placement$y[string_placement$gene == "RAD51C"] <- .665


# Subset to genes in network
# tmp.mtx.rad.blm <- as.matrix(gene.mtx[rownames(gene.mtx) %in% string_placement$gene,
#                                       colnames(gene.mtx) %in% string_placement$gene])

tmp.mtx.rad.blm <- as.matrix(gene.mtx[rownames(gene.mtx) %in% rad51.fil,
                                      colnames(gene.mtx) %in% rad51.fil])

mtx.melt <- melt(tmp.mtx.rad.blm)
mtx.melt$wt <- 0
mtx.melt$wt[mtx.melt$value > 1.5274] <- 1
mtx.melt$wt[mtx.melt$value < -1.5274] <- -1
mtx.melt$wt[mtx.melt$value < -2.754] <- -2
mtx.melt$wt[mtx.melt$value > 2.7555] <- 2
mtx.melt$First <- apply(mtx.melt[,1:2], 1, min)
mtx.melt$Second <- apply(mtx.melt[,1:2], 1, max)
mtx.melt$withinclust <- !((mtx.melt$First %in% blm & mtx.melt$Second %in% rad51.fil) |
                            (mtx.melt$First %in% rad51.fil & mtx.melt$Second %in% blm))
mtx.melt$inrad <- mtx.melt$First %in% rad51.fil & mtx.melt$Second %in% rad51.fil
mtx.melt$inblm <- mtx.melt$First %in% blm & mtx.melt$Second %in% blm
mtx.melt2 <- unique(mtx.melt[mtx.melt$wt != 0,c("First", "Second", "value", "wt", "withinclust", "inrad", "inblm")])

# Adjust ones acting weird in graph
mtx.melt2$First[mtx.melt2$First == "NUDC" & mtx.melt2$Second == "PSMC3IP"] <- "PSMC3IP"
mtx.melt2$Second[mtx.melt2$First == "PSMC3IP" & mtx.melt2$Second == "PSMC3IP"] <- "NUDC"
mtx.melt2$First[mtx.melt2$First == "RAD51D" & mtx.melt2$Second == "SWI5"] <- "SWI5"
mtx.melt2$Second[mtx.melt2$First == "SWI5" & mtx.melt2$Second == "SWI5"] <- "RAD51D"
mtx.melt2$First[mtx.melt2$First == "BRCA2" & mtx.melt2$Second == "SWI5"] <- "SWI5"
mtx.melt2$Second[mtx.melt2$First == "SWI5" & mtx.melt2$Second == "SWI5"] <- "BRCA2"

# Get information together for string map
string_interactions_short$X.node1[string_interactions_short$X.node1 == "NSMCE3"] <- "NDNL2"
string_interactions_short$int <- paste(string_interactions_short$X.node1, string_interactions_short$node2, sep = ":")
string_interactions_short$int2 <- paste(string_interactions_short$node1_string_id, string_interactions_short$node2_string_id, sep = ":")
string_link2 <- string_links[string_links$protein1 %in% string_interactions_short$node1_string_id & 
                               string_links$protein2 %in% string_interactions_short$node2_string_id,]
string_link2$int <- paste(string_link2$protein1, string_link2$protein2, sep = ":")
string_link3 <- string_link2[string_link2$int %in% string_interactions_short$int2[string_interactions_short$combined_score > .4],]
string_link3 <- inner_join(string_link3, string_interactions_short[,c("int", "int2")], by = c("int" = "int2"))

string_interactions_short$humanscreen <- string_interactions_short$int %in% string_link3$int.y[string_link3$experiments > 0][string_link3$int.y[string_link3$experiments > 0] %in% mtx.melt2$int]

string_mtx <- string_interactions_short[string_interactions_short$combined_score >= .4,c("X.node1", "node2", "combined_score", "humanscreen")]
string_mtx$wt <- ifelse(string_mtx$combined_score >= .9, 2, 1)
string_mtx$inrad <- string_mtx$X.node1 %in% rad51.fil & string_mtx$node2 %in% rad51.fil
string_mtx$inblm <- string_mtx$X.node1 %in% blm & string_mtx$node2 %in% blm
string_mtx$screen <- string_mtx$wt
string_mtx$screen[string_mtx$humanscreen] <- 3

string_mtx2 <- string_mtx[string_mtx$X.node1 %in% rad51.fil & string_mtx$node2 %in% rad51.fil,]

# Make graphs for the two datasets
gr <- graph_from_data_frame(mtx.melt2, directed = TRUE,
                            vertices = string_placement$gene[string_placement$gene %in% rad51.fil]) %>%
  as_tbl_graph() %>%
  mutate(N = centrality_degree(mode = 'in')) 
  
#gr <- as_tbl_graph(mtx.melt2) %>% 
# mutate(N = centrality_degree(mode = 'in')) 
gr <- gr %>%
  set_vertex_attr("type", value = vertex_attr(gr, "name") %in% blm)

gr.string <- graph_from_data_frame(string_mtx2, directed = TRUE, 
                                   vertices = string_placement$gene[string_placement$gene %in% rad51.fil]) %>%
  as_tbl_graph() %>%
  mutate(N = centrality_degree(mode = 'in')) 

#gr.string <- as_tbl_graph(string_mtx) %>% 
#  mutate(N = centrality_degree(mode = 'in')) 
gr.string <- gr.string %>%
  set_vertex_attr("type", value = vertex_attr(gr.string, "name") %in% blm)

# Make df for node placements
l1 <- data.frame(n = string_placement$gene[string_placement$gene %in% rad51.fil],
                 x = string_placement$x[string_placement$gene %in% rad51.fil],
                 y = string_placement$y[string_placement$gene %in% rad51.fil])
l1$n <- factor(l1$n, levels = V(gr)$name)
l1 <- l1[order(l1$n),]

l2 <- data.frame(n = string_placement$gene[string_placement$gene %in% rad51.fil],
                 x = string_placement$x[string_placement$gene %in% rad51.fil],
                 y = string_placement$y[string_placement$gene %in% rad51.fil])

# l2 <- data.frame(n = string_placement$gene[!string_placement$gene %in% c("NUDC", "WDR48")],
#                  x = string_placement$x[!string_placement$gene %in% c("NUDC", "WDR48")],
#                  y = string_placement$y[!string_placement$gene %in% c("NUDC", "WDR48")])
l2$n <- factor(l2$n, levels = V(gr.string)$name)
l2 <- l2[order(l2$n),]


g1 <- gr %>%
  ggraph(layout = l1) + 
  geom_edge_link(aes(width = as.character(abs(wt)), col = wt > 0)) + 
  geom_node_point(col = "#a7a9ac", size = 9) + 
  geom_node_text(aes(label = name), 
                 repel = FALSE, col = "black",
                 size = 2, fontface = "bold") +
  scale_edge_width_manual(values = c(.45, 1.25)) +
  scale_edge_color_manual(values = c("dodgerblue", "darkorange1")) +
  scale_color_viridis(discrete = TRUE) +
  scale_edge_alpha_manual(values = c(0, 1)) +
  theme_graph() +
  theme(legend.position = "none")

g2 <- gr.string %>%
  ggraph(layout = l2) + 
  geom_edge_link(aes(col = as.character(wt), width = as.character(abs(wt)))) + 
  geom_node_point(col = "#a7a9ac", size = 9) + 
  geom_node_text(aes(label = name), 
                 repel = FALSE, col = "black",
                 size = 2, fontface = "bold") +
  scale_edge_width_manual(values = c(.45, 1.25)) +
  scale_color_viridis(discrete = TRUE) +
  scale_edge_alpha_manual(values = c(0, 1)) +
  scale_edge_color_manual(values = c("#DDAA33", "#009988", "#BB5566")) +
  theme_graph() +
  theme(legend.position = "none")

pdf("../FIGURES/FIGURE-3/rad51_blm_interaction_network_stringnodeplacement.pdf", height = 5, width = 5)
g1
dev.off()

pdf("../FIGURES/FIGURE-3/rad51_string_network_stringnodeplacement.pdf", height = 5, width = 5)
g2
dev.off()
