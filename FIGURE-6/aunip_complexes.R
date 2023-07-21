library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("helper_functions.R")

con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_avg.txt", sep = ""))
con.tau.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_tau_oi_avg.txt", sep = ""))
gene.id.map <- read_tsv("../DATA/genecombination_id_map.txt")

rho.single.pheno <- read_tsv("../DATA/single_sgRNA_rho_phenotypes.txt")
rho.single.pheno$gene <- gsub("_.*", "", rho.single.pheno$sgRNA.ID)

rho.single.gene.pheno <- rho.single.pheno %>%
  group_by(gene) %>%
  summarise(rho = mean(Rho.Avg))

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

gi.cmp2$Highlight <- grepl("^AUNIP:|:AUNIP$", gi.cmp2$GeneCombinationName)
gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]

gi.cmp2$Complex <- NA
gi.cmp2$Complex[gi.cmp2$Highlight & gi.cmp2$Gene1 %in% c("ATR", "ATRIP")] <- "ATR-ATRIP complex"
gi.cmp2$Complex[gi.cmp2$Highlight & gi.cmp2$Gene2 %in% c("RAD9A", "RAD1", "HUS1")] <- "9-1-1 complex"
gi.cmp2$Complex[gi.cmp2$Highlight & gi.cmp2$Gene2 %in% c("FANCA", "FANCB", "FANCC", "FANCE",
                                                         "FANCF", "FANCG", "FANCL", "FANCM",
                                                         "C17orf70", "C19orf40", "STRA13")] <- "Fanconi anemia nuclear complex"
gi.cmp2$Complex[gi.cmp2$Highlight & gi.cmp2$Gene2 %in% c("BABAM1", "BRCA1", "BRCC3", "FAM175A", "UIMC1")] <- "BRCA1-A complex"

gi.cmp2$Non <- NA
gi.cmp2$Non[gi.cmp2$Highlight] <- gi.cmp2$Gene2[gi.cmp2$Highlight]
gi.cmp2$Non[gi.cmp2$Highlight & gi.cmp2$Non == "AUNIP"] <- gi.cmp2$Gene1[gi.cmp2$Highlight & 
                                                                                   gi.cmp2$Non == "AUNIP"]

gi.cmp3 <- left_join(gi.cmp2, rho.single.gene.pheno,
                     by = c("Non" = "gene"))

aunip.nolab <- ggplot(gi.cmp3, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.5, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.5, linetype= "dashed") +
  geom_point(alpha = .5, size = .85,  col = "#999999", shape = 16) +
  geom_point(data = gi.cmp3[gi.cmp3$Highlight & gi.cmp3$Non != "PARP1",],#aes(fill = rho), 
             size = 2, alpha = .5, col = "black", shape = 21) +
  geom_point(data = gi.cmp3[gi.cmp3$Highlight & gi.cmp3$Non == "PARP1",],col = "black", size = 2.2) +
  geom_point(data = gi.cmp3[!is.na(gi.cmp3$Complex),],
             aes(color = Complex), size = 3, shape = 15) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) + 
  coord_fixed() +
  scale_x_continuous(limits = c(-6, 4.5), breaks = seq(-6, 3, 3)) +
  scale_y_continuous(limits = c(-12,12), breaks = seq(-12, 12, 3)) +
  scale_fill_gradientn(colors = c("#003C30", "#01665E",  "#35978F", "#80CDC1", 
                                  "#E6E6E6", "#8C510A"),
                       values = c(0, .16, .32, .48, .66, .81, 1),
                       breaks = seq(-.9, .3, .3),
                       labels = c(-.9, -.6, -.3, 0, .3)) +
  scale_color_manual(values = c("#DDAA33", "#EE3377", "#ff2500", "#33bbee")) +
  xlab("") + ylab("")

aunip.lab <- ggplot(gi.cmp3[gi.cmp3$Highlight,], aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
  geom_abline(alpha = .7) +
  geom_abline(slope = 1, intercept = nu.bound, alpha = 0.5, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_hline(yintercept = tau.bound, alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -1*tau.bound, alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = gamma.bound, alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -1*(gamma.bound), alpha= 0.5, linetype= "dashed") +
  geom_point(alpha = .5, size = .85,  col = "#999999", shape = 16) +
  geom_point(data = gi.cmp3[gi.cmp3$Highlight & gi.cmp3$Non != "PARP1",],aes(fill = rho), 
             size = 2, alpha = .5, col = "black", shape = 21) +
  geom_point(data = gi.cmp3[gi.cmp3$Highlight & gi.cmp3$Non == "PARP1",],col = "black", size = 2.2) +
  geom_point(data = gi.cmp3[!is.na(gi.cmp3$Complex),],
             aes(color = Complex), size = 2.2, shape = 15) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) + 
  coord_fixed() +
  scale_x_continuous(limits = c(-6, 4.5), breaks = seq(-6, 3, 3)) +
  scale_y_continuous(limits = c(-12,12), breaks = seq(-12, 12, 3)) +
  scale_fill_gradientn(colors = c("#003C30", "#01665E",  "#35978F", "#80CDC1", 
                                  "#E6E6E6", "#8C510A"),
                       values = c(0, .16, .32, .48, .66, .81, 1),
                       breaks = seq(-.9, .3, .3),
                       labels = c(-.9, -.6, -.3, 0, .3)) +
  scale_color_manual(values = c("#DDAA33", "#EE3377", "#ff2500", "#33bbee")) +
  annotate(geom = "text", x = -6, y = 9, hjust = 0,
           label = paste0("r = ", round(cor(gi.cmp3$InteractionScore.Tau[gi.cmp3$Highlight],
                                            gi.cmp3$rho[gi.cmp3$Highlight]), 3))) +
  xlab("Gamma") + ylab("Tau")

ggsave("../FIGURES/FIGURE-6/aunip_complexes.png", aunip.nolab,
       device = "png", height = 10, width = 6, dpi = 300)
ggsave("../FIGURES/FIGURE-6/aunip_complexes_labels.svg", aunip.lab,
       device = "svg", height = 10, width = 6)
