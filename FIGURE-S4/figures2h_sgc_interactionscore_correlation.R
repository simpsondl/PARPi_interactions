library(ggplot2)
library(ggExtra)
library(readr)
library(dplyr)

# Load data
# Construct interaction scores
con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.r1 <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_r1.txt", sep = ""))
con.gamma.r2 <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_r2.txt", sep = ""))

# Merge replicate datasets
cols <- c(colnames(con.gamma.r1)[1:10], "GI.z")
con.gamma <- full_join(con.gamma.r1[,colnames(con.gamma.r1) %in% cols], 
                       con.gamma.r2[,colnames(con.gamma.r2) %in% cols],
                       by = colnames(con.gamma.r1)[1:10],
                       suffix = c(".R1", ".R2"))
con.gamma.sgc <- con.gamma %>% group_by(GuideCombinationID) %>%
  summarise(GI.z.R1 = mean(GI.z.R1), GI.z.R2 = mean(GI.z.R2))

con.gamma.sgc <- inner_join(con.gamma.sgc, unique(con.gamma[,c("GuideCombinationID", "Category")]))

# Rearrange data for plotting
tmp.df <- rbind(con.gamma.sgc[con.gamma.sgc$Category == "X+Y",], 
                con.gamma.sgc[con.gamma.sgc$Category == "X+NT",],
                con.gamma.sgc[con.gamma.sgc$Category == "NT+NT",],
                con.gamma.sgc[con.gamma.sgc$Category == "X+X",])

# Sample subset for plotting labeled pngs
tmp.df2 <- tmp.df[sample(nrow(tmp.df),10000),]
# Check that the sampling chose at least one member of each category
stopifnot(length(table(tmp.df2$Category)) == 4)

# Label-free theme
nolabel_theme <-   theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.title.y = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.grid = element_blank(),
                         axis.ticks.length = unit(.2, "cm"),
                         legend.position = "none")
nolabel_axes <- list(xlab(""), ylab(""), ggtitle("")) 

# No labels for flat PNG - FIGURE S2H
ps2h <- ggplot(tmp.df,aes(GI.z.R1, GI.z.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .35, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  theme_bw() + 
  nolabel_theme +
  nolabel_axes +
  scale_x_continuous(limit = c(-18,13.5), breaks = c(-15, -10, -5, 0, 5, 10)) +
  scale_y_continuous(limit = c(-18,13.5), breaks = c(-15, -10, -5, 0, 5, 10))

# Add marginal histograms - FIGURE S2H
ps2hh <- ggMarginal(ps2h, groupColour = FALSE, groupFill = TRUE)

# labels for SVG - FIGURE S2H
ps2hl <- ggplot(tmp.df2,aes(GI.z.R1, GI.z.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .5, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  annotate("text", hjust = 0, x= -10, y = 6, 
           label = paste("r =", round(cor(tmp.df$GI.z.R1[tmp.df$Category %in% c("X+X", "X+Y")], tmp.df$GI.z.R2[tmp.df$Category %in% c("X+X", "X+Y")]),3))) +
  theme_bw() +
  removeGrid() +
  xlab("Gamma Interaction Score Replicate 1") +
  ylab("Gamma Interaction Score Replicate 2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_continuous(limit = c(-18,13.5), breaks = c(-15, -10, -5, 0, 5, 10)) +
  scale_y_continuous(limit = c(-18,13.5), breaks = c(-15, -10, -5,0, 5, 10))

# Save plots
ggsave("../FIGURES/FIGURE-S4/gamma_sgRNA_score_replicate_correlation.png", 
       ps2hh, device = "png", width = 2.5, height = 2.5, dpi = 300)
ggsave("../FIGURES/FIGURE-S4/gamma_sgRNA_score_replicate_correlation_labels.svg", 
       ps2hl, device = "svg", width = 10, height = 10, dpi = 300)
