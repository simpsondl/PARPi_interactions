library(readr)
library(ggplot2)
library(dplyr)

# Load functions
source("helper_functions.R")

# Load data
gene_path_pfx <- "../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/"
string <- read_tsv("../DATA/External_Data/string_interactions.tsv")
colnames(string)[1] <- "node1"
gamma.gene <- read_tsv(paste0(gene_path_pfx,
                              "gene_combination_interaction_scores_gamma_oi_avg.txt"))
id.map <- read_tsv("../DATA/genecombination_id_map.txt")

# Update string results to match IDs in score tables
string$node1[string$node1 == "BHLHE40"] <- "STRA13"
string$node2[string$node2 == "BHLHE40"] <- "STRA13"
# Set up GeneCombinationName column
string$First <- apply(string[,1:2], 1, min)
string$Second <- apply(string[,1:2], 1, max)
string$GeneCombinationName <- paste(string$First, string$Second, sep = ":")
# Merge in gc ids
string_merge <- inner_join(string[,c("GeneCombinationName","combined_score")],
                           id.map) %>% unique()
gamma.gene <- inner_join(gamma.gene, id.map)
# Identify validated interactions
gamma.gene$Validated <- gamma.gene$GeneCombinationID %in% string_merge$GeneCombinationID

# Remove all but nontargeting combinations
gamma.gene <- gamma.gene[gamma.gene$Category %in% c("X+X", "X+Y"),]

# Make a 538 x 538 grid with gene names
gene.grid <- expand.grid(unique(gamma.gene$Gene1), unique(gamma.gene$Gene1))
# Merge in interaction scores
gene.grid <- left_join(gene.grid, 
                       gamma.gene, 
                       by = c("Var1" = "Gene2", "Var2" = "Gene1"))
gene.grid <- left_join(gene.grid, 
                       gamma.gene, 
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
# Rename rows
rownames(gene.mtx) <- gene.mtx$Var1
# Remove redundant column
gene.mtx <- gene.mtx[,2:ncol(gene.mtx)]

gene.mtx <- impute.knn(as.matrix(gene.mtx))$data

# Calculate correlations
gene.cors <- cor(as.matrix(gene.mtx))

# Set up for correlation plots
cor.melt <- melt(gene.cors)
cor.melt$First <- apply(cor.melt[,1:2], 1, min)
cor.melt$Second <- apply(cor.melt[,1:2], 1, max)
cor.melt$GeneCombinationName <- paste(cor.melt$First, cor.melt$Second, sep = ":")
cor.melt2 <- inner_join(cor.melt, id.map[,c("GeneCombinationName", "GeneCombinationID")] )
cor.melt2 <- cor.melt2[cor.melt2$First != cor.melt2$Second,]
cor.melt2$Validated <- cor.melt2$GeneCombinationID %in% string_merge$GeneCombinationID
cor.melt2$HighConfidence <- cor.melt2$GeneCombinationID %in% string_merge$GeneCombinationID[string_merge$combined_score >= .7]
cor.melt2$HighestConfidence <- cor.melt2$GeneCombinationID %in% string_merge$GeneCombinationID[string_merge$combined_score >= .9]
cor.melt3 <- unique(cor.melt2[,c(7:10,3)])
cor.melt3$Confidence[cor.melt3$HighestConfidence] <- "Highest"
cor.melt3$Confidence[!cor.melt3$HighestConfidence & (cor.melt3$Validated)] <- "Medium/High"
cor.melt3$Confidence[is.na(cor.melt3$Confidence)] <- "No"
cor.melt3$Confidence <- factor(cor.melt3$Confidence, levels = rev(c("Highest", "Medium/High", "No")))
cor.melt3$Significance <- "NS"
cor.melt3$Significance[cor.melt3$GeneCombinationID %in% 
                     gamma.gene$GeneCombinationID[abs(gamma.gene$InteractionScore) > 2.723]] <- "Strongly Significant"

cor.melt3$Highlight <- as.character(cor.melt3$Confidence)
tmp <- cor.melt3[cor.melt3$Significance != "NS",]
tmp$Highlight <- "Strongly Significant"

cor.melt3 <- rbind(cor.melt3, tmp)

cor.melt3$Highlight <- factor(cor.melt3$Highlight, levels = c("No", "Medium/High", "Highest", "Strongly Significant"))
cor.melt3 <- cor.melt3[order(cor.melt3$Confidence),]

cor.plt <- ggplot(cor.melt3, aes(value, col = Highlight)) + 
  geom_density(size = 1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Gamma Gene-Gene Profile Correlation") +
  ylab("Density") +
  scale_color_manual(values = c("#999999", "#DDAA33", "#009988", "#BB5566")) +
  annotate(geom = "text", x = .2, y = 3.5, hjust = 0,
           label = "Unknown:") +
  annotate(geom = "text", x = .2, y = 3.32, hjust = 0,
           label = "Medium/High Confidence:") +
  annotate(geom = "text", x = .2, y = 3.14, hjust = 0,
           label = "Highest Confidence:") +
  annotate(geom = "text", x = .2, y = 2.96, hjust = 0,
           label = "Strongly Significant:") +
  annotate(geom = "text", x = .73, y = 3.5, hjust = 1,
           label = table(cor.melt3$Highlight)[[1]]) +
  annotate(geom = "text", x = .73, y = 3.32, hjust = 1,
           label = table(cor.melt3$Highlight)[[2]]) +
  annotate(geom = "text", x = .73, y = 3.14, hjust = 1,
           label = table(cor.melt3$Highlight)[[3]]) +
  annotate(geom = "text", x = .73, y = 2.96, hjust = 1,
           label = table(cor.melt3$Highlight)[[4]]) 

ggsave("../FIGURES/FIGURE-2/figure_2e_stringdb_correlation_confidences_density.svg",
       cor.plt, device = "svg", width = 7, height = 4, dpi = 300)

