library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(impute)

# Load functions
source("helper_functions.R")

# Load data
gc_path_pfx <- "../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/"
gamma.r1 <- read_tsv(paste0(gc_path_pfx,
                                  "gene_combination_interaction_scores_gamma_oi_r1.txt"))
gamma.r2 <- read_tsv(paste0(gc_path_pfx,
                            "gene_combination_interaction_scores_gamma_oi_r2.txt"))
gamma.gene <- read_tsv(paste0(gc_path_pfx,
                              "gene_combination_interaction_scores_gamma_oi_avg.txt"))
#string <- read_tsv("../DATA/External_Data/string_interactions.tsv")
#colnames(string)[1] <- "node1"
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
clusts <- read_tsv("../DATA/heatmap_clusters_power4power6.txt")

fancs <- clusts$gene[clusts$P6 == 2]

# Update string results to match IDs in score tables
#string$node1[string$node1 == "BHLHE40"] <- "STRA13"
#string$node2[string$node2 == "BHLHE40"] <- "STRA13"
# Set up GeneCombinationName column
#string$First <- apply(string[,1:2], 1, min)
#string$Second <- apply(string[,1:2], 1, max)
#string$GeneCombinationName <- paste(string$First, string$Second, sep = ":")
# Merge in gc ids
#string_merge <- inner_join(string[,c("GeneCombinationName","combined_score")],
#                           id.map) %>% unique()

# gamma.gene <- inner_join(gamma.gene, id.map)
# # Identify validated interactions
# #gamma.gene$Validated <- gamma.gene$GeneCombinationID %in% string_merge$GeneCombinationID
# 
# # Remove all but X+Y combinations
# gamma.gene <- gamma.gene[gamma.gene$Category %in% c("X+Y"),]
# 
# # Make a 538 x 538 grid with gene names
# gene.grid <- expand.grid(unique(gamma.gene$Gene1), unique(gamma.gene$Gene1))
# # Merge in interaction scores
# gene.grid <- left_join(gene.grid, 
#                        gamma.gene, 
#                        by = c("Var1" = "Gene2", "Var2" = "Gene1"))
# gene.grid <- left_join(gene.grid, 
#                        gamma.gene, 
#                        by = c("Var1" = "Gene1", "Var2" = "Gene2"))
# # Combine the columns
# gene.grid$Gene.GI <- ifelse(is.na(gene.grid$InteractionScore.x), gene.grid$InteractionScore.y, gene.grid$InteractionScore.x)
# gene.grid$N <- ifelse(is.na(gene.grid$N.x), gene.grid$N.y, gene.grid$N.x)
# 
# # Remove extraneous columns
# gene.grid <- gene.grid[,c("Var1", "Var2", "Gene.GI", "N")]
# gene.grid$Gene.GI <- as.numeric(gene.grid$Gene.GI)
# gene.grid$N <- as.numeric(gene.grid$N)
# 
# #gene.grid.imputed <- impute_missing_genegene(gene.grid)
# 
# # Go from long to wide
# gene.mtx <- pivot_wider(gene.grid[,1:3], names_from = Var2, values_from = Gene.GI)
# gene.mtx <- as.data.frame(gene.mtx)
# ### Rename rows
# rownames(gene.mtx) <- gene.mtx$Var1
# ### Remove redundant column
# gene.mtx <- gene.mtx[,2:ncol(gene.mtx)]
# 
# # Impute missing values
# gene.mtx <- impute.knn(as.matrix(gene.mtx))$data
# 
# # Calculate correlations
#gene.cors <- cor(as.matrix(gene.mtx))

# Set up for correlation plots
# cor.melt <- melt(gene.cors)
# cor.melt$First <- apply(cor.melt[,1:2], 1, min)
# cor.melt$Second <- apply(cor.melt[,1:2], 1, max)
# cor.melt$GeneCombinationName <- paste(cor.melt$First, cor.melt$Second, sep = ":")
# cor.melt2 <- inner_join(cor.melt, id.map[,c("GeneCombinationName", "GeneCombinationID")] )
# cor.melt2 <- cor.melt2[cor.melt2$First != cor.melt2$Second,]
# cor.melt2$Validated <- cor.melt2$GeneCombinationID %in% string_merge$GeneCombinationID
# cor.melt2$HighConfidence <- cor.melt2$GeneCombinationID %in% string_merge$GeneCombinationID[string_merge$combined_score >= .7]
# cor.melt2$HighestConfidence <- cor.melt2$GeneCombinationID %in% string_merge$GeneCombinationID[string_merge$combined_score >= .9]
# cor.melt3 <- unique(cor.melt2[,c(7:10,3)])
# cor.melt3$Confidence[cor.melt3$HighestConfidence] <- "Highest"
# cor.melt3$Confidence[is.na(cor.melt3$Confidence)] <- "No"
# cor.melt3$Confidence <- factor(cor.melt3$Confidence, levels = rev(c("Highest", "No")))
# cor.melt3 <- cor.melt3[order(cor.melt3$Confidence),]
# cor.melt4 <- left_join(cor.melt3, id.map)

# Merge in gene names
gamma.r1 <- inner_join(gamma.r1,
                       id.map)
gamma.r2 <- inner_join(gamma.r2,
                       id.map)
gamma.gene <- inner_join(gamma.gene,
                         id.map)

gamma.gene$fanc <- gamma.gene$Gene1 %in% fancs & gamma.gene$Gene2 %in% fancs

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol.gamma.r1 <- gamma.r1[!grepl("non-targeting", gamma.r1$GeneCombinationName),]
noncontrol.gamma.r2 <- gamma.r2[!grepl("non-targeting", gamma.r2$GeneCombinationName),]

# Make a 538 x 538 grid with gene names
gene.grid.r1 <- expand.grid(unique(noncontrol.gamma.r1$Gene1), unique(noncontrol.gamma.r1$Gene1))
gene.grid.r2 <- expand.grid(unique(noncontrol.gamma.r2$Gene1), unique(noncontrol.gamma.r2$Gene1))
# Merge in interaction scores
gene.grid.r1 <- left_join(gene.grid.r1, 
                       noncontrol.gamma.r1, 
                       by = c("Var1" = "Gene2", "Var2" = "Gene1"))
gene.grid.r1 <- left_join(gene.grid.r1, 
                       noncontrol.gamma.r1, 
                       by = c("Var1" = "Gene1", "Var2" = "Gene2"))
gene.grid.r2 <- left_join(gene.grid.r2, 
                          noncontrol.gamma.r2, 
                          by = c("Var1" = "Gene2", "Var2" = "Gene1"))
gene.grid.r2 <- left_join(gene.grid.r2, 
                          noncontrol.gamma.r2, 
                          by = c("Var1" = "Gene1", "Var2" = "Gene2"))
# Combine the columns
gene.grid.r1$Gene.GI <- ifelse(is.na(gene.grid.r1$InteractionScore.x), gene.grid.r1$InteractionScore.y, gene.grid.r1$InteractionScore.x)
gene.grid.r1$N <- ifelse(is.na(gene.grid.r1$N.x), gene.grid.r1$N.y, gene.grid.r1$N.x)
gene.grid.r2$Gene.GI <- ifelse(is.na(gene.grid.r2$InteractionScore.x), gene.grid.r2$InteractionScore.y, gene.grid.r2$InteractionScore.x)
gene.grid.r2$N <- ifelse(is.na(gene.grid.r2$N.x), gene.grid.r2$N.y, gene.grid.r2$N.x)

# Remove extraneous columns
gene.grid.r1 <- gene.grid.r1[,c("Var1", "Var2", "Gene.GI", "N")]
gene.grid.r1$Gene.GI <- as.numeric(gene.grid.r1$Gene.GI)
gene.grid.r1$N <- as.numeric(gene.grid.r1$N)
gene.grid.r2 <- gene.grid.r2[,c("Var1", "Var2", "Gene.GI", "N")]
gene.grid.r2$Gene.GI <- as.numeric(gene.grid.r2$Gene.GI)
gene.grid.r2$N <- as.numeric(gene.grid.r2$N)

#gene.grid.r1.imputed <- impute_missing_genegene(gene.grid.r1)
#gene.grid.r2.imputed <- impute_missing_genegene(gene.grid.r2)

# Go from long to wide
gene.mtx.r1 <- pivot_wider(gene.grid.r1[,1:3], names_from = Var2, values_from = Gene.GI)
gene.mtx.r1 <- as.data.frame(gene.mtx.r1)
gene.mtx.r2 <- pivot_wider(gene.grid.r2[,1:3], names_from = Var2, values_from = Gene.GI)
gene.mtx.r2 <- as.data.frame(gene.mtx.r2)
# Rename rows
rownames(gene.mtx.r1) <- gene.mtx.r1$Var1
rownames(gene.mtx.r2) <- gene.mtx.r2$Var1
# Remove redundant column
gene.mtx.r1 <- gene.mtx.r1[,2:ncol(gene.mtx.r1)]
gene.mtx.r2 <- gene.mtx.r2[,2:ncol(gene.mtx.r2)]

# Impute missing values
gene.mtx.r1 <- impute.knn(as.matrix(gene.mtx.r1))$data
gene.mtx.r2 <- impute.knn(as.matrix(gene.mtx.r2))$data

# Calculate correlations
gene.cor.r1 <- cor(as.matrix(gene.mtx.r1))
gene.cor.r2 <- cor(as.matrix(gene.mtx.r2))
# Wide to long
cors.r1 <- melt(gene.cor.r1)
cors.r2 <- melt(gene.cor.r2)
cors <- inner_join(cors.r1, cors.r2,
                   by = c("Var1", "Var2"),
                   suffix = c(".R1",".R2"))
cors2 <- cors[cors$Var1 != cors$Var2,]
# Remove duplicates
cors2$First <- apply(cors2[,c(1,2)], 1, min)
cors2$Second <- apply(cors2[,c(1,2)], 1, max)
cors2$GeneCombinationName <- paste(cors2$First, cors2$Second, sep = ":")
cors3 <- inner_join(cors2, id.map)
cors3$fanc <- cors3$Gene1 %in% fancs & cors3$Gene2 %in% fancs
cors4 <- unique(cors3[,c("GeneCombinationID", "value.R1", "value.R2", "fanc")])

# Merge in string confidence level
#cors4 <- left_join(cors4, cor.melt3[,c("GeneCombinationID", "Confidence")])

# Add interaction significance
cors4$Significance <- "NS"
cors4$Significance[cors4$GeneCombinationID %in% 
                     gamma.gene$GeneCombinationID[gamma.gene$InteractionScore >= 2.7555 | 
                                                    gamma.gene$InteractionScore <= -2.754]] <- "Strongly Significant"
cors4$Significance[cors4$fanc] <- "FA pathway"
# Combine interaction and string confidence for label
# cors4$Highlight <- as.character(cors4$Significance)
# cors4$Highlight[cors4$Confidence != "No"] <- as.character(cors4$Confidence[cors4$Confidence != "No"])
# cors4$Highlight[cors4$Highlight == "Significant"] <- "Other"
# cors4$Highlight[cors4$Highlight == "Strongly Significant"] <- "|Gamma| >= 2.723"
# cors4$Highlight[cors4$Highlight == "Highest"] <- "String >= .9"


# Rearrange for plotting
cors.tmp <- rbind(cors4[cors4$Significance == "NS",],
                  cors4[cors4$Significance == "Strongly Significant",],
                  cors4[cors4$Significance == "FA pathway",])

cors.tmp$Significance <- factor(cors.tmp$Significance, levels = c("NS", "Strongly Significant", "FA pathway"))

# cors.tmp <- rbind(cors4[cors4$Highlight == "Other",],
#                   cors4[cors4$Highlight == "|Gamma| >= 2.723",],
#                   cors4[cors4$Highlight == "String >= .9",])

cors.samp <- cors4[sample(1:nrow(cors4), 10000),]

# Label-free theme
nolabel_theme <-   theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.title.y = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.grid = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_line(),
                         axis.ticks.length = unit(.2, "cm"),
                         legend.position = "none")
nolabel_axes <- list(xlab(""), ylab(""), ggtitle("")) 

# Plot
# No labels for flat PNG
plt <- ggplot(cors.tmp, aes(value.R1, value.R2, col = Significance)) +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .7, size = .25, pch = 16) +
  geom_point(data = cors.tmp[cors.tmp$Significance == "Strongly Significant",], alpha = .7, size = .3, pch = 16) +
  geom_point(data = cors.tmp[cors.tmp$Significance == "FA pathway",], alpha = .7, size = .3, pch = 16) +
  theme_bw() + 
  nolabel_axes +
  nolabel_theme +
  scale_color_manual(values = c("#999999", "#BB5566", "#33bbee")) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-.2, .6, .2)) + 
  scale_y_continuous(breaks = seq(-.2, .6, .2))

plt.lab <- ggplot(cors.samp, aes(value.R1, value.R2, col = Significance)) +
        geom_hline(yintercept = 0, alpha = .5) +
        geom_vline(xintercept = 0, alpha = .5) +
        geom_point(alpha = .3, size = 1, pch = 16) +
        theme_bw() + 
        removeGrid() +
        xlab("Replicate 1") + ylab("Replicate 2") +
        scale_color_manual(values = c("#888888", "#BB5566", "#33bbee")) +
        coord_fixed() +
        scale_x_continuous(breaks = seq(-.2, .6, .2)) + 
        scale_y_continuous(breaks = seq(-.2, .6, .2)) + 
        guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
        annotate(geom = "text", x = .3, y = -.25, label = paste0("r = ", 
                                                                 round(cor(cors4$value.R1, cors4$value.R2), 3))) +
        annotate(geom = "text", x = .3, y = -.3, label = paste0("r = ", 
                                                                 round(cor(cors4$value.R1[cors4$Significance == "Strongly Significant"], 
                                                                           cors4$value.R2[cors4$Significance == "Strongly Significant"]), 3)))

# Save plots
ggsave("../FIGURES/FIGURE-2/gene_gi_replicate_correlation_scatterplot.png", plt, 
       device = "png", width = 2.5, height = 2.5, dpi = 300)
ggsave("../FIGURES/FIGURE-2/figure_2d_gene_gi_replicate_correlation_scatterplot_labels.svg", plt.lab, 
       device = "svg", width = 8, height = 8, dpi = 300)
