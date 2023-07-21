phenos <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
rhos <- read_tsv("../DATA/single_sgRNA_rho_phenotypes.txt")
g <- read_tsv("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/gene_combination_interaction_scores_gamma_oi_avg.txt")
t <- read_tsv("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/gene_combination_interaction_scores_tau_oi_avg.txt")
n <- read_tsv("../DATA/Interaction_Scores/Compiled/GeneCombination_Scores/gene_combination_interaction_scores_nu_oi_avg.txt")
idmap <- read_tsv("../DATA/genecombination_id_map.txt")

phenos$gene <- gsub("_.*", "", phenos$sgRNA.ID)
rhos$gene <- gsub("_.*", "", rhos$sgRNA.ID)
pgenes <- phenos %>% group_by(gene) %>% summarise(Gamma = mean(Gamma.Avg), Tau = mean(Tau.Avg))
rgenes <- rhos %>% group_by(gene) %>% summarise(Rho = mean(Rho.Avg))

g <- inner_join(g, idmap)
t <- inner_join(t, idmap)
n <- inner_join(n, idmap)

g.genes <- setdiff(unique(c(g$Gene1, g$Gene2)), "non-targeting")
t.genes <- setdiff(unique(c(t$Gene1, t$Gene2)), "non-targeting")
n.genes <- intersect(g.genes, t.genes)

df <- data.frame(gene = unique(c(g.genes, t.genes, n.genes)),
                 ng.tot = 0, nt.tot = 0, nn.tot = 0)
for(i in 1:nrow(df)){
  df$ng.tot[i] <- sum(abs(g$InteractionScore[g$Gene1 == df$gene[i] | g$Gene2 == df$gene[i]]) >= 1.5274)
  df$nt.tot[i] <- sum(abs(t$InteractionScore[t$Gene1 == df$gene[i] | t$Gene2 == df$gene[i]]) >= 1.5301)
  df$nn.tot[i] <- sum(abs(n$InteractionScore[n$Gene1 == df$gene[i] | n$Gene2 == df$gene[i]]) >= 2.977)
}

df <- left_join(df, pgenes)
df <- left_join(df, rgenes)

gplt <- ggplot(df[df$gene %in% g.genes & df$Gamma <= -.05,], aes(Gamma, ng.tot)) +
  #geom_vline(xintercept = -.05, linetype = "dotted", alpha = .5) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  scale_x_continuous(breaks = seq(-.1, -.3, -.1)) +
  ylab("Number of significant interactions") +
  annotate(geom = "text", x = -.35, y = 175, hjust = 0,
           label = paste0("r = ", round(cor(df$Gamma[df$gene %in% g.genes & df$Gamma <= -.05],
                                            df$ng.tot[df$gene %in% g.genes & df$Gamma <= -.05]), 3)))

tplt <- ggplot(df[df$gene %in% t.genes & df$Tau <= -.05,], aes(Tau, nt.tot)) +
  #geom_vline(xintercept = -.05, linetype = "dotted", alpha = .5) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  scale_x_continuous(breaks = seq(-.2, -.6, -.2)) +
  ylab("Number of significant interactions") +
  annotate(geom = "text", x = -.6, y = 175, hjust = 0,
           label = paste0("r = ", round(cor(df$Tau[df$gene %in% t.genes & df$Tau <= -.05],
                                            df$nt.tot[df$gene %in% t.genes & df$Tau <= -.05]), 3)))

nplt <- ggplot(df[df$gene %in% n.genes & abs(df$Rho) >= .05 & df$gene != "PARP1",], aes(abs(Rho), nn.tot)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  #scale_x_continuous(breaks = seq(-.1, -.7, -.1),
  #                   labels = c("", -0.2, "", -0.4, "", -0.6, "")) +
  ylab("Number of significant interactions") +
  annotate(geom = "text", x = .1, y = 55, hjust = 0,
           label = paste0("r = ", round(cor(abs(df$Rho[df$gene %in% n.genes & abs(df$Rho) >= .05 & df$gene != "PARP1"]),
                                            df$nn.tot[df$gene %in% n.genes & abs(df$Rho) >= .05 & df$gene != "PARP1"]), 3)))


ggsave("../FIGURES/FIGURE-S4/gammaphenotype_v_nsigints_scatterplot.svg", gplt,
       device = "svg", height = 4, width = 4)
ggsave("../FIGURES/FIGURE-S4/tauphenotype_v_nsigints_scatterplot.svg", tplt,
       device = "svg", height = 4, width = 4)
ggsave("../FIGURES/FIGURE-S4/rhophenotype_v_nsigints_scatterplot.svg", nplt,
       device = "svg", height = 4, width = 4)
