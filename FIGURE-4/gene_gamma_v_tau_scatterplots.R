library(readr)
library(dplyr)
library(ggplot2)

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

gi.cmp2$Highlight <- grepl("^PARP1:|:PARP1$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

parp1 <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8,  col = "#999999") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,],col = "black", size = 1.2) +
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
  xlab("") + ylab("") 

ggsave("../FIGURES/FIGURE-4/parp1_highlight.png", parp1, width = 6, height = 10, dpi = 300)

gi.cmp2$Highlight <- grepl("^RIF1:|:RIF1$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

rif1 <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8,  col = "#999999") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,],col = "black", size = 1.2) +
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
  xlab("") + ylab("") 

ggsave("../FIGURES/FIGURE-4/rif1_highlight.png", rif1, width = 6, height = 10, dpi = 300)

gi.cmp2$Highlight <- grepl("^FANCA:|:FANCA$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

fanca <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8,  col = "#999999") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,],col = "black", size = 1.2) +
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
  xlab("") + ylab("") 

ggsave("../FIGURES/FIGURE-4/fanca_highlight.png", fanca, width = 6, height = 10, dpi = 300)

gi.cmp2$Highlight <- grepl("^DNPH1:|:DNPH1$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

dnph1 <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8,  col = "#999999") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,],col = "black", size = 1.2) +
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
  xlab("") + ylab("") 

ggsave("../FIGURES/FIGURE-4/dnph1_highlight.png", dnph1, width = 6, height = 10, dpi = 300)

gi.cmp2$Highlight <- grepl("^BRCA2:|:BRCA2$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

brca2 <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8,  col = "#999999") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,],col = "black", size = 1.2) +
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
  xlab("") + ylab("") 

ggsave("../FIGURES/FIGURE-4/brca2_highlight.png", brca2, width = 6, height = 10, dpi = 300)

gi.cmp2$Highlight <- grepl("^AUNIP:|:AUNIP$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

aunip <-ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8,  col = "#999999") +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight,],col = "black", size = 1.2) +
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
  xlab("") + ylab("") 

ggsave("../FIGURES/FIGURE-4/aunip_highlight.png", aunip, width = 6, height = 10, dpi = 300)

