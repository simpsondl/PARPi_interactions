library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)

# Load data
phenotypes <- read_tsv("../DATA/raw_phenotypes.txt")

brca2_combos <- phenotypes$GeneCombinationID[phenotypes$FirstGene == "BRCA2" |
                                               phenotypes$SecondGene == "BRCA2"]

outliers <- c("BARD1", "BRCA1", "BRCA2", "DCTN5", "HUS1", 
              "LRRC45", "MED23", "PALB2", "PARP1", "RAD17", 
              "RAD54L", "RAD9A", "SEPSECS", "XRCC1")
outlier_combos <- phenotypes$GeneCombinationID[phenotypes$FirstGene %in% outliers | 
                                                 phenotypes$SecondGene %in% outliers]

# Generate gene level phenotypes
pheno <- phenotypes %>% 
  group_by(GeneCombinationID) %>% 
  summarise(Gamma.R1 = mean(Gamma.R1), Gamma.R2 = mean(Gamma.R2))

pheno$BRCA2 <- pheno$GeneCombinationID %in% brca2_combos
pheno$Outlier <- pheno$GeneCombinationID %in% outlier_combos

# Reorder for plotting
pheno2 <- rbind(pheno[!(pheno$BRCA2 & pheno$Outlier),],
                pheno[pheno$Outlier,],
                pheno[pheno$BRCA2,])

# Downsample for label output
tmp.df <- pheno2[sample(nrow(pheno2),10000),]
tmp.df$combined <- "Other"
tmp.df$combined[tmp.df$Outlier] <- "Outlier"
tmp.df$combined[tmp.df$BRCA2] <- "BRCA2"

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

plt <- ggplot(pheno2, aes(Gamma.R1, Gamma.R2)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .6, alpha = .3, pch = 16, col = "#bbbbbb") +
  geom_point(data = pheno2[pheno2$Outlier,], col = "#888888", size= .6, alpha = .3, pch = 16) +
  geom_point(data = pheno2[pheno2$BRCA2,], size = .8, alpha = .7) +
  xlim(c(-.6,.15)) +
  ylim(c(-.6,.15)) +
  theme_bw() + 
  nolabel_axes + 
  nolabel_theme

plt.lab <- ggplot(tmp.df, aes(Gamma.R1, Gamma.R2, col = combined)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .5, alpha = .3, pch = 16) +
  scale_color_manual(values = c("black", "#bbbbbb", "#888888")) +
  theme_bw() + 
  removeGrid() + 
  xlab("Replicate 1") + 
  ylab(" Replicate 2") + 
  xlim(c(-.6,.15)) +
  ylim(c(-.6,.15)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  annotate("text", x = -.5, y = .08, label = paste("r =", round(cor(pheno2$Gamma.R1, pheno2$Gamma.R2), 3)))

ggsave("../FIGURES/TechnicalNote/combination_phenotype_outliers_highlighted.png", plt,
       device = "png", height = 6, width = 6, dpi = 300)
ggsave("../FIGURES/TechnicalNote/combination_phenotype_outliers_highlighted_labels.pdf", plt.lab,
       device = "pdf", height = 6, width = 6, dpi = 300)
