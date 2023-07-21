library(readr)
library(ggplot2)

# Load functions
source("helper_functions.R")

con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_avg.txt", sep = ""))
con.tau.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_tau_oi_avg.txt", sep = ""))
gene.id.map <- read_tsv("../DATA/genecombination_id_map.txt")

con.oi.avg <- compute_construct_differential_scores(con.gamma.oi.avg, con.tau.oi.avg)

gamma.gene.gi <- compute_gene_interaction_scores(con.gamma.oi.avg, "GI.z")
tau.gene.gi <- compute_gene_interaction_scores(con.tau.oi.avg, "GI.z")
nu.gene.gi <- compute_gene_interaction_scores(con.oi.avg, "DGI.z")

gamma.bound <- 1.5274
gamma.hc.pos.bound <- 2.7555
gamma.hc.neg.bound <- -2.754
tau.bound <- 1.5301
tau.hc.bound <- 2.883
nu.bound <- 2.977

gi.cmp <-inner_join(gamma.gene.gi[,1:3], tau.gene.gi[,1:3], by = c("GeneCombinationID", "Category"), suffix = c(".Gamma", ".Tau"))
gi.cmp <-inner_join(gi.cmp, nu.gene.gi[,1:3], by = c("GeneCombinationID", "Category"))
colnames(gi.cmp)[5] <- "InteractionScore.Nu"

gi.cmp2 <- inner_join(gi.cmp, gene.id.map, 
                      by = c("GeneCombinationID")) %>%
  filter(Category == "X+Y")

gi.cmp2$Highlight <- NA
gi.cmp2$Highlight[gi.cmp2$InteractionScore.Nu >= nu.bound] <- "Resistance Interaction"
gi.cmp2$Highlight[gi.cmp2$InteractionScore.Nu <= -1*nu.bound] <- "Sensitizing Interaction"
gi.cmp2$Highlight[is.na(gi.cmp2$Highlight)] <- "No Interaction"
gi.cmp2$Highlight <- factor(gi.cmp2$Highlight, levels = c("No Interaction",
                                                          "Sensitizing Interaction", "Resistance Interaction"))

gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]
gi.samp <- gi.cmp2[sample(1:nrow(gi.cmp2),10000),]

plt <- ggplot(gi.cmp2[gi.cmp2$Category == "X+Y",], aes(InteractionScore.Gamma, InteractionScore.Tau)) +
  geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8) +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight != "No Interaction",], aes(color = Highlight), alpha = .5, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "None") + coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("") + ylab("") +
  scale_color_manual(values = c("#33716B", "#D81B60" ))

plt_lab <- ggplot(gi.samp[gi.samp$Category == "X+Y",], aes(InteractionScore.Gamma, InteractionScore.Tau)) +
  geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
  geom_point(alpha = .5, size = .8) +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight != "No Interaction",], aes(color = Highlight), alpha = .5, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) + coord_fixed() +
  scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
  scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
  xlab("Gamma IS") + ylab("Tau IS") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  scale_color_manual(values = c("#33716B", "#D81B60" ))

ggsave("../FIGURES/FIGURE-4/gamma_tau_interaction_category_scatterplot.png", plt,
       device = "png", height = 10, width = 6)
ggsave("../FIGURES/FIGURE-4/gamma_tau_interaction_category_scatterplot_labels.svg", plt_lab,
       device = "svg", height = 10, width = 6)
