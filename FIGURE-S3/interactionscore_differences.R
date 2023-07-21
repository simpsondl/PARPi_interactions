library(readr)
library(dplyr)
library(ggplot2)
library(ggpattern)

gamma.r1.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_r1.txt")
gamma.r2.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_r2.txt")
gamma.avg.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")

id.map <- read_tsv("../DATA/genecombination_id_map.txt")

# Load data
phenos <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")

# Extract gene name from sgRNA.ID
phenos$gene <- gsub("_.*", "", phenos$sgRNA.ID)
phenos2 <- phenos %>% group_by(gene) %>% summarise(mean = mean(Gamma.Avg))


outliers <- c("BARD1", "BRCA1", "BRCA2", "DCTN5", "HUS1", 
              "LRRC45", "MED23", "PALB2", "PARP1", "RAD17", 
              "RAD54L", "RAD9A", "SEPSECS", "XRCC1")
double.outlier <- c("BARD1", "BRCA2", "MED23", "PARP1", "RAD54L")

all.gamma <- inner_join(gamma.r1.all[,c(1:10,13,16)], 
                        gamma.r2.all[,c(1:10,13,16)], 
                        by = colnames(gamma.r1.all)[1:10], 
                        suffix = c(".R1", ".R2"))
all.gamma <- inner_join(all.gamma, gamma.avg.all[,c(1:10,13,16)], 
                        by = colnames(gamma.r1.all)[1:10])

outlier.ints <- id.map$GeneCombinationID[id.map$Gene1 %in% outliers | id.map$Gene2 %in% outliers]
double.outlier.ints <- id.map$GeneCombinationID[id.map$Gene1 %in% double.outlier | id.map$Gene2 %in% double.outlier]
outlier.gamma <- all.gamma[all.gamma$GeneCombinationID %in% outlier.ints,]
double.outlier.gamma <- all.gamma[all.gamma$GeneCombinationID %in% double.outlier.ints,]

double.outlier.gamma$BARD <- double.outlier.gamma$GeneCombinationID %in% 
                              id.map$GeneCombinationID[id.map$Gene1 == "BARD1" | id.map$Gene2 == "BARD1"]
double.outlier.gamma$BRCA <- double.outlier.gamma$GeneCombinationID %in% 
  id.map$GeneCombinationID[id.map$Gene1 == "BRCA2" | id.map$Gene2 == "BRCA2"]
double.outlier.gamma$MED <- double.outlier.gamma$GeneCombinationID %in% 
  id.map$GeneCombinationID[id.map$Gene1 == "MED23" | id.map$Gene2 == "MED23"]
double.outlier.gamma$PARP <- double.outlier.gamma$GeneCombinationID %in% 
  id.map$GeneCombinationID[id.map$Gene1 == "PARP1" | id.map$Gene2 == "PARP1"]
double.outlier.gamma$RAD <- double.outlier.gamma$GeneCombinationID %in% 
  id.map$GeneCombinationID[id.map$Gene1 == "RAD54L" | id.map$Gene2 == "RAD54L"]

double.int.df <- double.outlier.gamma[double.outlier.gamma$BARD,c(4,9:16)]
double.int.df$Gene <- "BARD1"

tmp <- double.outlier.gamma[double.outlier.gamma$BRCA, c(4,9:16)]
tmp$Gene <- "BRCA2"
double.int.df <- rbind(double.int.df, tmp)

tmp <- double.outlier.gamma[double.outlier.gamma$MED, c(4,9:16)]
tmp$Gene <- "MED23"
double.int.df <- rbind(double.int.df, tmp)

tmp <- double.outlier.gamma[double.outlier.gamma$PARP, c(4,9:16)]
tmp$Gene <- "PARP1"
double.int.df <- rbind(double.int.df, tmp)

tmp <- double.outlier.gamma[double.outlier.gamma$RAD, c(4,9:16)]
tmp$Gene <- "RAD54L"
double.int.df <- rbind(double.int.df, tmp)

double.m <- melt(double.int.df[,c(5,7,9,10)])

bx.plt <- ggplot(double.m, aes(Gene, value)) + 
  geom_violin_pattern(aes(pattern = variable), pattern_spacing = .015) +
  geom_boxplot(aes(group = interaction(Gene, variable)), 
               outlier.shape = NA, coef = 0,
               position = position_dodge(.9), width = .2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  ylab("Interaction Score") +
  scale_pattern_discrete(name = "Replicate", labels = c("R1", "R2", "AVG"))

ggsave("../FIGURES/TechnicalNote/outliergene_interactionscore_distributions.pdf", bx.plt,
       device = "pdf", height = 6, width = 6, dpi = 300)

gene.cors <- data.frame(Gene = unique(id.map$Gene1), Cor = NA)
for(i in gene.cors$Gene){
  tmp <- all.gamma[all.gamma$GeneCombinationID %in% id.map$GeneCombinationID[id.map$Gene1 == i | id.map$Gene2 == i],]
  if(nrow(tmp) > 0){
    gene.cors$Cor[gene.cors$Gene == i] <- cor(tmp$GI.z.R1, tmp$GI.z.R2, use = "complete.obs")
  } 
}
gene.cors <- gene.cors[!is.na(gene.cors$Cor) & gene.cors$Gene != "non-targeting",]
gene.cors$Highlight <- gene.cors$Gene %in% double.outlier

genecor.plt <- ggplot(gene.cors, aes(Cor)) + 
  geom_vline(xintercept = mean(gene.cors$Cor), linetype = "dashed", alpha = .5) +
  geom_density(fill = "#bbbbbb", alpha = .5) + 
  geom_text(data = gene.cors[gene.cors$Highlight,], aes(y = .3, label = Gene), size = 3) + 
  geom_rug(alpha = .5) +
  geom_rug(data = gene.cors[gene.cors$Highlight,], col = "red", length = unit(.05, "npc"), size = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  ylab("") + xlab("Replicate Correlation")

ggsave("../FIGURES/TechnicalNote/gene_replicatecorrelation.pdf", genecor.plt,
       device = "pdf", height = 6, width = 6, dpi = 300)

gamma.gcs <- all.gamma %>% group_by(GeneCombinationID) %>% 
  summarise(GI.R1 = mean(GI.z.R1), GI.R2 = mean(GI.z.R2), GI.Avg = mean(GI.z), N = n())
gamma.gcs <- inner_join(gamma.gcs, id.map)

sig.df <- data.frame(gene = unique(gamma.gcs$Gene1), n = 0)
for(i in unique(gamma.gcs$Gene1)){
  sig.df$n[sig.df$gene == i] <- sum(abs(gamma.gcs$GI.Avg[gamma.gcs$Gene1 == i | gamma.gcs$Gene2 == i]) > 1.5274)
}

cor.sig <- inner_join(gene.cors, sig.df, by = c("Gene" = "gene"))
cor.fit <- lm(n ~ Cor, data = cor.sig)
cor.sig$diff <- cor.sig$n - (cor.fit$coefficients[[2]]*cor.sig$Cor + cor.fit$coefficients[[1]])

outlier.cors <- rosnerTest(cor.sig$diff, k = 20)$all.stats
cor.sig$Highlight <- cor.sig$Gene %in% cor.sig$Gene[outlier.cors$Obs.Num[outlier.cors$Outlier]]

cor.sig.plt <- ggplot(cor.sig) + 
  theme_bw() + 
  geom_abline(slope = cor.fit$coefficients[[2]],
              intercept = cor.fit$coefficients[[1]]) +
  geom_abline(slope = cor.fit$coefficients[[2]],
              intercept = cor.fit$coefficients[[1]] + 57.5, linetype = "dashed") +
  geom_abline(slope = cor.fit$coefficients[[2]],
              intercept = cor.fit$coefficients[[1]] - 57.5, linetype = "dashed") +
  geom_point(aes(Cor, n, size = Highlight, alpha = Highlight)) + 
  geom_point(data = cor.sig[cor.sig$Gene %in% outliers,],aes(Cor, n), col =  "#DDAA33", size = 2) +
  geom_point(data = cor.sig[cor.sig$Gene %in% double.outlier,],aes(Cor, n), col =  "#DDAA33", size = 2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  xlab("Replicate Correlation") + ylab("Number of Significant Interactions") +
  annotate(geom = "text", x = .5, y = 15, hjust = 0,
           label = paste("r =",round(cor(cor.sig$n, cor.sig$Cor), 3))) +
  scale_alpha_manual(values = c(.7, 1)) +
  scale_size_manual(values = c(1.3,2))

ggsave("../FIGURES/TechnicalNote/numsiginteractions_vs_replicatecorrelation.pdf", cor.sig.plt,
       device = "pdf", height = 6, width = 6, dpi = 300)
ggsave("../FIGURES/TechnicalNote/numsiginteractions_vs_replicatecorrelation.svg", cor.sig.plt,
       device = "svg", height = 6, width = 6)
