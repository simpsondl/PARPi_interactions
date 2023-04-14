library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)

# Load functions
source("helper_functions.R")

con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.r1 <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_r1.txt", sep = ""))
con.gamma.r2 <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_r2.txt", sep = ""))
gene.id.map <- read_tsv("../DATA/genecombination_id_map.txt")
pgene.id.map <- read_tsv("../DATA/pseudogenecombination_id_map.txt")
gcs <- read_tsv("../DATA/nulldist_pvals_gamma_genegeneinteractionscores.txt")

gamma.gene.gi.r1 <- compute_pseudogene_interaction_scores(con.gamma.r1, "GI.z")
gamma.gene.gi.r2 <- compute_pseudogene_interaction_scores(con.gamma.r2, "GI.z")

# Merge replicate datasets
gc.gamma <- full_join(gamma.gene.gi.r1[,c("PseudogeneCombinationID", "Category", "InteractionScore")], 
                      gamma.gene.gi.r2[,c("PseudogeneCombinationID", "Category", "InteractionScore")],
                      by = c("PseudogeneCombinationID", "Category"),
                      suffix = c(".R1", ".R2"))

# Map gene combination ids to pseudogenecombination ids
id.map <- inner_join(gene.id.map[,3:4], pgene.id.map[,3:4],
                     by = c("GeneCombinationName" = "PseudogeneCombinationName"))
#id.map$Sig <- id.map$GeneCombinationID %in% gcs$X1[gcs$Qvalue <= .05]
id.map$Sig <- id.map$GeneCombinationID %in% gcs$X1[gcs$LocalFDR <= .01]

gc.gamma$Sig <- gc.gamma$PseudogeneCombinationID %in% id.map$PseudogeneCombinationID[id.map$Sig]

gc.gamma$CatSig <- gc.gamma$Category
gc.gamma$CatSig[gc.gamma$Sig & gc.gamma$Category != "X+X"] <- gc.gamma$Sig[gc.gamma$Sig & gc.gamma$Category != "X+X"]

# Rearrange data for plotting
tmp.df <- rbind(gc.gamma[gc.gamma$CatSig == "X+Y",],
                gc.gamma[gc.gamma$CatSig == "X+NT",],
                gc.gamma[gc.gamma$CatSig == "NT+NT",],
                gc.gamma[gc.gamma$CatSig == "TRUE",])

# tmp.df <- rbind(gc.gamma[gc.gamma$Category == "X+Y",],
#                 gc.gamma[gc.gamma$Category == "X+NT",],
#                 gc.gamma[gc.gamma$Category == "NT+NT",])

tmp.df$CatSig <- factor(tmp.df$CatSig, levels = c("X+Y", "X+NT", "NT+NT", "TRUE"))
tmp.df$Category <- factor(tmp.df$Category, levels = c("X+Y", "X+NT", "NT+NT"))

# Sample subset for plotting labeled pngs
tmp.df2 <- tmp.df[sample(nrow(tmp.df),10000),]
tmp.df2 <- rbind(tmp.df2,
                 tmp.df[tmp.df$InteractionScore.R1 == min(tmp.df$InteractionScore.R1),])
tmp.df2 <- rbind(tmp.df2,
                 tmp.df[tmp.df$InteractionScore.R2 == min(tmp.df$InteractionScore.R2),])
tmp.df2 <- rbind(tmp.df2,
                 tmp.df[tmp.df$InteractionScore.R1 == max(tmp.df$InteractionScore.R1),])
tmp.df2 <- rbind(tmp.df2,
                 tmp.df[tmp.df$InteractionScore.R2 == max(tmp.df$InteractionScore.R2),])
# Check that the sampling chose both categories
stopifnot(length(table(tmp.df2$CatSig)) == 4)

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

# No labels for flat PNG - FIGURE 2C
p2c <- ggplot(tmp.df,aes(InteractionScore.R1, InteractionScore.R2, color = CatSig)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_hline(yintercept = -1.524, alpha = .5, linetype = "dashed") +
  geom_hline(yintercept = 1.524, alpha = .5, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_vline(xintercept = -1.524, alpha = .5, linetype = "dashed") +
  geom_vline(xintercept = 1.524, alpha = .5, linetype = "dashed") +
  geom_point(size= 1, alpha = .6, pch = 16) +
  #geom_point(data = tmp.df[tmp.df$Sig,], color = "#BB5566", size= 1, alpha = .2, pch = 16) +
  scale_color_manual(values = c("#999999", "#abddde", "#046c9a", "#BB5566")) +
  theme_bw() + 
  nolabel_theme +
  nolabel_axes +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-9, 9, 3)) +
  scale_y_continuous(breaks = seq(-12, 6, 3))

# Add marginal histograms - FIGURE 2C
p2ch <- ggMarginal(p2c, groupColour = FALSE, groupFill = TRUE)

# labels for SVG - FIGURE 2C
p2cl <- ggplot(tmp.df2,aes(InteractionScore.R1, InteractionScore.R2, color = CatSig)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= 1, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "#999999", "#BB5566")) +
  annotate("text", hjust = 0, x= -8, y = 6, 
           label = paste("r[a] ==", round(cor(tmp.df$InteractionScore.R1, 
                                          tmp.df$InteractionScore.R2),3)), parse = TRUE) +
  annotate("text", hjust = 0, x= -8, y = 5.2, 
           label = paste("r[s] ==", round(cor(tmp.df$InteractionScore.R1[tmp.df$CatSig %in% c("TRUE")], 
                                          tmp.df$InteractionScore.R2[tmp.df$CatSig %in% c("TRUE")]),3)), parse = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Gamma Gene-Gene Interaction Score Replicate 1") +
  ylab("Gamma Gene-Gene Interaction Score Replicate 2") +
  labs(color = "") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-9, 9, 3)) +
  scale_y_continuous(breaks = seq(-12, 6, 3))

# Save plots
ggsave("../FIGURES/FIGURE-2/figure_2c_genegene_gammainteractionscore_replicate_correlation.png", 
       p2ch, device = "png", width = 8, height = 8, dpi = 300)
ggsave("../FIGURES/FIGURE-2/figure_2c_genegene_gammainteractionscore_replicate_correlation_labels.svg", 
       p2cl, device = "svg", width = 8, height = 8, dpi = 300)

