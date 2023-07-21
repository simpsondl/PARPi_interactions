library(ggplot2)
library(ggExtra)
library(readr)

# Load functions
source("helper_functions.R")

# Load GI screen data
phenos <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
rho <- read_tsv("../DATA/single_sgRNA_rho_phenotypes.txt")
phenos <- inner_join(phenos, rho)
# Extract gene name from sgRNA.ID
phenos$Gene <- gsub("_.*", "", phenos$sgRNA.ID)
# Group by gene
gene.phenos <- phenos %>% 
                group_by(Gene) %>% 
                summarise(Gamma.R1 = mean(Gamma.R1), Gamma.R2 = mean(Gamma.R2), Gamma.Avg = mean(Gamma.Avg),
                          Tau.R1 = mean(Tau.R1), Tau.R2 = mean(Tau.R2), Tau.Avg = mean(Tau.Avg),
                          Rho.R1 = mean(Rho.R1), Rho.R2 = mean(Rho.R2), Rho.Avg = mean(Rho.Avg))

# Load primary screen data and calculate phenotypes
primary_counts <- read_delim("../DATA/raw_counts_primary_screen.txt", 
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE, trim_ws = TRUE, skip = 4)

colnames(primary_counts) <- c("sgRNA.ID","NIRAP.R1","NIRAP.R2","T0.R1","T0.R2", "DMSO.R1", "DMSO.R2")

# Define conditions in each arm
cond.df <- data.frame(Colname = colnames(primary_counts[,2:7]),
                      Samplename = gsub("\\.R.*", "", colnames(primary_counts[,2:7])),
                      Replicate = gsub(".*\\.","", colnames(primary_counts[,2:7])))

# Calculate phenotypes for primary screen
raw_phenotypes <- calculate_primary_screen_phenotypes(counts = primary_counts, 
                                       conds = cond.df,
                                       pseudocount = 10,
                                       normalize = TRUE,
                                       doublings = c(7.65, 7.76, 4.71, 6.07))
# Extract gene name
raw_phenotypes$Gene <- gsub("_.*", "", raw_phenotypes$sgRNA.ID)
# Group by Gene
ps_gene_gamma <- raw_phenotypes %>% 
            group_by(Gene) %>% 
            slice_max(order_by = abs(Gamma.Avg), n = 3) %>%
            summarise(Gamma.R1 = mean(Gamma.R1), Gamma.R2 = mean(Gamma.R2), Gamma.Avg = mean(Gamma.Avg))

ps_gene_tau <- raw_phenotypes %>% 
  group_by(Gene) %>% 
  slice_max(order_by = abs(Tau.Avg), n = 3) %>%
  summarise(Tau.R1 = mean(Tau.R1), Tau.R2 = mean(Tau.R2), Tau.Avg = mean(Tau.Avg))

ps_gene_rho <- raw_phenotypes %>% 
  group_by(Gene) %>% 
  slice_max(order_by = abs(Rho.Avg), n = 3) %>%
  summarise(Rho.R1 = mean(Rho.R1), Rho.R2 = mean(Rho.R2), Rho.Avg = mean(Rho.Avg))

ps_gene_phenos <- inner_join(ps_gene_gamma, ps_gene_tau)
ps_gene_phenos <- inner_join(ps_gene_phenos, ps_gene_rho)

# Merge primary and interaction screen phenotypes
both_screens <- inner_join(gene.phenos, ps_gene_phenos,
                           by = "Gene",
                           suffix = c(".GI",".P"))
# Fit model
fitg <- lm(Gamma.Avg.GI ~ Gamma.Avg.P, data = both_screens)
fitt <- lm(Tau.Avg.GI ~ Tau.Avg.P, data = both_screens)
fitr <- lm(Rho.Avg.GI ~ Rho.Avg.P, data = both_screens)

# Make plots
plt.g <- ggplot(both_screens, aes(Gamma.Avg.P, Gamma.Avg.GI)) + 
  geom_point(size=.8, alpha = .7) +
  theme_bw() + 
  removeGrid() +
  geom_abline(slope = fitg$coefficients[[2]], intercept = fitg$coefficients[[1]], col = "red") + 
  # annotate(geom = "text", x = -.4, y = .05, hjust = 0,
  #          label = paste("y = ", round(fitg$coefficients[[2]],3), "x + ", round(fitg$coefficients[[1]],3), sep = "")) +
  annotate(geom = "text", x = -.1, y = -.35, hjust = 0,
           label = paste0("r = ", round(cor(both_screens$Gamma.Avg.GI, both_screens$Gamma.Avg.P), 3))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Gene Gamma Phenotype - Primary Screen") + 
  ylab("Gene Gamma Phenotype - Interaction Screen") 

plt.t <- ggplot(both_screens, aes(Tau.Avg.P, Tau.Avg.GI)) + 
  geom_point(size=.8, alpha = .7) +
  theme_bw() + 
  removeGrid() +
  geom_abline(slope = fitt$coefficients[[2]], intercept = fitt$coefficients[[1]], col = "red") + 
  annotate(geom = "text", x = -.5, y = .15, hjust = 0,
           label = paste("y = ", round(fitt$coefficients[[2]],3), "x + ", round(fitt$coefficients[[1]],3), sep = "")) +
  annotate(geom = "text", x = .1, y = -.62, hjust = 0,
           label = paste0("r = ", round(cor(both_screens$Tau.Avg.GI, both_screens$Tau.Avg.P), 3))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Gene Tau Phenotype - Primary Screen") + 
  ylab("Gene Tau Phenotype - Interaction Screen") 

plt.r <- ggplot(both_screens, aes(Rho.Avg.P, Rho.Avg.GI)) + 
  geom_point(size=.8, alpha = .7) +
  theme_bw() + 
  removeGrid() +
  geom_abline(slope = fitr$coefficients[[2]], intercept = fitr$coefficients[[1]], col = "red") + 
  annotate(geom = "text", x = -.5, y = .15, hjust = 0,
           label = paste("y = ", round(fitr$coefficients[[2]],3), "x + ", round(fitr$coefficients[[1]],3), sep = "")) +
  annotate(geom = "text", x = .1, y = -.62, hjust = 0,
           label = paste0("r = ", round(cor(both_screens$Rho.Avg.GI, both_screens$Rho.Avg.P), 3))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Gene Rho Phenotype - Primary Screen") + 
  ylab("Gene Rho Phenotype - Interaction Screen") 


# Save plots
ggsave("../FIGURES/FIGURE-S2/screen_comparison_gene_gamma_correlation_scatterplot.svg", 
       plt.g, device = "svg", width = 5, height = 5)
ggsave("../FIGURES/FIGURE-S2/screen_comparison_gene_tau_correlation_scatterplot.svg", 
       plt.t, device = "svg", width = 5, height = 5)
ggsave("../FIGURES/FIGURE-S2/screen_comparison_gene_rho_correlation_scatterplot.svg", 
       plt.r, device = "svg", width = 5, height = 5)
