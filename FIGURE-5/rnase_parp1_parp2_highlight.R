library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("helper_functions.R")

con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_avg.txt", sep = ""))
con.tau.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_tau_oi_avg.txt", sep = ""))
gene.id.map <- read_tsv("../DATA/genecombination_id_map.txt")

con.oi.avg <- compute_construct_differential_scores(con.gamma.oi.avg, con.tau.oi.avg)

gamma.gene.gi <- compute_gene_interaction_scores(con.gamma.oi.avg, "GI.z")
tau.gene.gi <- compute_gene_interaction_scores(con.tau.oi.avg, "GI.z")
nu.gene.gi <- compute_gene_interaction_scores(con.oi.avg, "DGI.z")

gamma.bound <- 1.524
gamma.hc.bound <- 2.723
tau.bound <- 1.527
tau.hc.bound <- 2.853
nu.bound <- 2.975

gi.cmp <-inner_join(gamma.gene.gi[,1:3], tau.gene.gi[,1:3], by = c("GeneCombinationID", "Category"), suffix = c(".Gamma", ".Tau"))
gi.cmp <-inner_join(gi.cmp, nu.gene.gi[,1:3], by = c("GeneCombinationID", "Category"))
colnames(gi.cmp)[5] <- "InteractionScore.Nu"

gi.cmp2 <- inner_join(gi.cmp, gene.id.map, 
                      by = c("GeneCombinationID")) %>%
  filter(Category == "X+Y")

genes <- c("RNASEH2A", "RNASEH2B", "RNASEH2C", "PARP1", "PARP2")
gis.to.highlight <- gene.id.map[gene.id.map$Gene1 %in% genes & gene.id.map$Gene2 %in% genes,]
gis.to.highlight <- gis.to.highlight[gis.to.highlight$Gene1 != gis.to.highlight$Gene2,]

gi.cmp2$Highlight <- gi.cmp2$GeneCombinationID %in% gis.to.highlight$GeneCombinationID
gi.cmp2$Col[gi.cmp2$Highlight & gi.cmp2$Gene1 == "PARP1" & grepl("RNASEH",gi.cmp2$Gene2)] <- "PARP1"
gi.cmp2$Col[gi.cmp2$Highlight & gi.cmp2$Gene1 == "PARP2" & grepl("RNASEH",gi.cmp2$Gene2)] <- "PARP2"
gi.cmp2$Col[gi.cmp2$Highlight & grepl("RNASEH",gi.cmp2$Gene1) & grepl("RNASEH",gi.cmp2$Gene2)] <- "RNASEH2A/B/C"

gi.cmp2.samp <- gi.cmp2[sample(1:nrow(gi.cmp2), 1000),]

rnaseh2.highlight <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .7,  col = "#bbbbbb") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,], aes(col = Col), size = 1.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none") + 
  coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("") + ylab("") +
  scale_color_manual(values = c( "#ddcc77", "#117733", "#CC6677"), 
                     na.value = "black", name = "",
                     breaks = c("PARP1", "PARP2", "RNASEH2A/B/C"))

rnaseh2.highlight.lab <-ggplot(gi.cmp2.samp, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .7,  col = "#bbbbbb") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,], aes(col = Col), size = 1.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) + 
  coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("Gamma") + ylab("Tau") +
  scale_color_manual(values = c("#500D50", "#999933", "#CC6677"),  
                     na.value = "black", name = "",
                     breaks = c("PARP1", "PARP2", "RNASEH2A/B/C"))
  
rnaseh2_names <- ggplot(gi.cmp2.samp, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .7,  col = "#bbbbbb") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,], aes(col = Col), size = 1.6) +
  geom_text_repel(data = gi.cmp2[gi.cmp2$Highlight,], aes(label = GeneCombinationName)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2, "cm")) + 
  coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("") + ylab("") +
  scale_color_manual(values = c("#500D50", "#999933", "#CC6677"), 
                     na.value = "black", name = "",
                     breaks = c("PARP1", "PARP2", "RNASEH2A/B/C"))

ggsave("../FIGURES/FIGURE-5/rnaseh2_parp1_parp2_highlight.png", rnaseh2.highlight, 
       device = "png", width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-5/rnaseh2_parp1_parp2_highlight_labels.svg", rnaseh2.highlight.lab, 
       device = "svg", width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-5/rnaseh2_parp1_parp2_highlight_labels_with_genenames.svg", rnaseh2_names, 
       device = "svg", width = 6, height = 10, dpi = 300)
