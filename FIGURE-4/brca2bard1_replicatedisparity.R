library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)

# Load data
phenotypes <- read_tsv("../DATA/raw_phenotypes.txt")

phenotypes$highlight <- phenotypes$FirstGene %in% c("BRCA2", "BARD1") | phenotypes$SecondGene %in% c("BRCA2", "BARD1")

# Generate gene level phenotypes
pheno <- phenotypes %>% 
  group_by(PseudogeneCombinationID) %>% 
  summarise(Gamma.R1 = mean(Gamma.R1), Gamma.R2 = mean(Gamma.R2), Category = unique(Category))


pheno$Highlight <- pheno$PseudogeneCombinationID %in% phenotypes$PseudogeneCombinationID[phenotypes$highlight]

# Reorder for plotting
pheno2 <- rbind(pheno[!pheno$Highlight,],
                pheno[pheno$Highlight,])

# Downsample for label output
tmp.df <- pheno2[sample(nrow(pheno2),10000),]
stopifnot(length(table(tmp.df$Category)) == 4)

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

plt <- ggplot(pheno2, aes(Gamma.R1, Gamma.R2, color = Highlight)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .9, alpha = .8, pch = 16) +
  #geom_point(data = pheno2[pheno2$Highlight,], size= 1.1, alpha = .7, pch = 16) +
  xlim(c(-.6,.15)) +
  ylim(c(-.6,.15)) +
  scale_color_manual(values = c("#999999", "black")) + 
  theme_bw() + 
  nolabel_axes + 
  nolabel_theme
plt.m <- ggMarginal(plt, groupColour = FALSE, groupFill = TRUE)
plt.lab <- ggplot(tmp.df, aes(Gamma.R1, Gamma.R2, color = Highlight)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .5, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#999999", "black")) + 
  theme_bw() + 
  removeGrid() + 
  xlab("Replicate 1") + 
  ylab(" Replicate 2") + 
  xlim(c(-.6,.15)) +
  ylim(c(-.6,.15)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  annotate("text", x = -.5, y = .08, label = paste("r =", round(cor(pheno2$Gamma.R1, pheno2$Gamma.R2), 3)))


ggsave("../FIGURES/FIGURE-4/figure_4_replicate_disparity_bard1brca2.png", plt.m, 
       device = "png", width = 8, height = 8, dpi = 300)
ggsave("../FIGURES/FIGURE-4/figure_4_replicate_disparity_bard1brca2_labels.svg", plt.lab, 
       device = "svg", width = 8, height = 8, dpi = 300)
