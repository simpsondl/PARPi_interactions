library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)

# Load data
phenotypes <- read_tsv("../DATA/raw_phenotypes.txt")
rho <- read_tsv("../DATA/raw_rho_phenotypes.txt")

phenotypes2 <- phenotypes %>% group_by(PseudogeneCombinationID) %>%
  summarise(Category = unique(Category),
            Gamma.R1 = mean(Gamma.R1),
            Gamma.R2 = mean(Gamma.R2),
            Tau.R1 = mean(Tau.R1),
            Tau.R2 = mean(Tau.R2))
rho2 <- rho %>% group_by(PseudogeneCombinationID) %>%
  summarise(Category = unique(Category),
            Rho.R1 = mean(Rho.R1),
            Rho.R2 = mean(Rho.R2))

phenotypes2 <- inner_join(phenotypes2,
                          rho2,
                          by = c("PseudogeneCombinationID", "Category"))

# Rearrange data for plotting
tmp.df <- rbind(phenotypes2[phenotypes2$Category == "X+Y",], 
                phenotypes2[phenotypes2$Category == "X+NT",],
                phenotypes2[phenotypes2$Category == "NT+NT",],
                phenotypes2[phenotypes2$Category == "X+X",])


# Sample subset for plotting labeled pngs
tmp.df2 <- tmp.df[sample(nrow(tmp.df),10000),]
# Check that the sampling chose at least one member of each category
stopifnot(length(table(tmp.df2$Category)) == 4)

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

# No labels for flat PNG - FIGURE S2C
g <- ggplot(tmp.df,aes(Gamma.R1, Gamma.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .4, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  scale_x_continuous(limits = c(-.6, .15), breaks = c(-.6, -.45, -.3, -.15, 0, .15)) +
  scale_y_continuous(limits = c(-.6, .15), breaks = c(-.6, -.45, -.3, -.15, 0, .15)) +
  theme_bw() + 
  nolabel_theme +
  nolabel_axes 

# Add marginal histograms - FIGURE S2C
gh <- ggMarginal(g, groupColour = FALSE, groupFill = TRUE)

# No labels for flat PNG - FIGURE S2C
t <- ggplot(tmp.df,aes(Tau.R1, Tau.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .4, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  scale_x_continuous(limits = c(-1, .35), breaks = c(-1, -.75, -.5, -.25, 0, .25)) +
  scale_y_continuous(limits = c(-1, .35), breaks = c(-1, -.75, -.5, -.25, 0, .25)) +
  theme_bw() + 
  nolabel_theme +
  nolabel_axes 

# Add marginal histograms - FIGURE S2C
th <- ggMarginal(t, groupColour = FALSE, groupFill = TRUE)

# No labels for flat PNG - FIGURE S2C
r <- ggplot(tmp.df,aes(Rho.R1, Rho.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .4, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  scale_x_continuous(limits = c(-1.85, 1.6), breaks = c(-1.5, -1, -.5, 0, .5, 1, 1.5)) +
  scale_y_continuous(limits = c(-1.85, 1.6), breaks = c(-1.5, -1, -.5, 0, .5, 1, 1.5)) +
  theme_bw() + 
  nolabel_theme +
  nolabel_axes 

# Add marginal histograms - FIGURE S2C
rh <- ggMarginal(r, groupColour = FALSE, groupFill = TRUE)

# labels for SVG - FIGURE 1A
gl <- ggplot(tmp.df2,aes(Gamma.R1, Gamma.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .5, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  annotate("text", hjust = 0, x= -.5, y = .1, 
           label = paste("r =", round(cor(tmp.df$Gamma.R1, tmp.df$Gamma.R2),3))) +
  scale_x_continuous(limits = c(-.6, .15), breaks = c(-.6, -.45, -.3, -.15, 0, .15)) +
  scale_y_continuous(limits = c(-.6, .15), breaks = c(-.6, -.45, -.3, -.15, 0, .15)) +
  theme_bw() +
  removeGrid() +
  xlab("Gamma Phenotype Replicate 1") +
  ylab("Gamma Phenotype Replicate 2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

tl <- ggplot(tmp.df2,aes(Tau.R1, Tau.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .5, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  annotate("text", hjust = 0, x= -.5, y = .1, 
           label = paste("r =", round(cor(tmp.df$Tau.R1, tmp.df$Tau.R2),3))) +
  scale_x_continuous(limits = c(-1.25, .4), breaks = c(-1.25, -1, -.75, -.5, -.25, 0, .25)) +
  scale_y_continuous(limits = c(-1.25, .4), breaks = c(-1.25, -1, -.75, -.5, -.25, 0, .25)) +
  theme_bw() +
  removeGrid() +
  xlab("Tau Phenotype Replicate 1") +
  ylab("Tau Phenotype Replicate 2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

rl <- ggplot(tmp.df2,aes(Rho.R1, Rho.R2, color = Category)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size= .5, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#046c9a", "#abddde", "darkorchid4", "#999999")) +
  annotate("text", hjust = 0, x= -.5, y = .1, 
           label = paste("r =", round(cor(tmp.df$Rho.R1, tmp.df$Rho.R2),3))) +
  scale_x_continuous(limits = c(-2.3, 2), breaks = c(-2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2)) +
  scale_y_continuous(limits = c(-2.3, 2), breaks = c(-2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2)) +
  theme_bw() +
  removeGrid() +
  xlab("Rho Phenotype Replicate 1") +
  ylab("Rho Phenotype Replicate 2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

# Save plots
ggsave("../FIGURES/FIGURE-S2/gammaphenotype_genereplicate_correlation.png", 
       gh, device = "png", width = 2.5, height = 2.5, dpi = 300)
ggsave("../FIGURES/FIGURE-S2/gammaphenotype_genereplicate_correlation_labels.svg", 
       gl, device = "svg", width = 6, height = 6)
ggsave("../FIGURES/FIGURE-S2/tauphenotype_genereplicate_correlation.png", 
       th, device = "png", width = 2.5, height = 2.5, dpi = 300)
ggsave("../FIGURES/FIGURE-S2/tauphenotype_genereplicate_correlation_labels.svg", 
       tl, device = "svg", width = 6, height = 6)
ggsave("../FIGURES/FIGURE-S2/rhophenotype_genereplicate_correlation.png", 
       rh, device = "png", width = 2.5, height = 2.5, dpi = 300)
ggsave("../FIGURES/FIGURE-S2/rhophenotype_genereplicate_correlation_labels.svg", 
       rl, device = "svg", width = 6, height = 6)
