library(ggplot2)
library(ggExtra)
library(readr)

# Load functions
source("helper_functions.R")

# Load annotations
annotations <- read_tsv("../DATA/External_Data/screen_annotations_w_mapid.txt")
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
ps_gene <- raw_phenotypes %>% 
  group_by(Gene) %>% 
  summarise(Gamma.R1 = mean(Gamma.R1), Gamma.R2 = mean(Gamma.R2), Gamma.Avg = mean(Gamma.Avg),
            Tau.R1 = mean(Tau.R1), Tau.R2 = mean(Tau.R2), Tau.Avg = mean(Tau.Avg))
ps_gene$Nu.Avg <- ps_gene$Tau.Avg - ps_gene$Gamma.Avg
ps_gene$Screen <- ps_gene$Gene %in% annotations$MapGeneName

plt.gamma <- ggplot(ps_gene[,c("Gene", "Gamma.Avg", "Screen")], aes(Gamma.Avg)) +
  geom_histogram(fill = "#999999", binwidth = .015) +
  geom_histogram(data = ps_gene[ps_gene$Screen,c("Gene", "Gamma.Avg", "Screen")], 
                 aes(Gamma.Avg), fill = "#9dcbe3", binwidth = .015) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Primary Screen Gamma Phenotypes") +
  ylab("Count") +
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(1, 10, 100, 1000, 10000))

plt.tau <- ggplot(ps_gene[,c("Gene", "Tau.Avg", "Screen")], aes(Tau.Avg)) +
  geom_histogram(fill = "#999999", binwidth = .02) +
  geom_histogram(data = ps_gene[ps_gene$Screen,c("Gene", "Tau.Avg", "Screen")], 
                 aes(Tau.Avg), fill = "#9dcbe3", binwidth = .02) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Primary Screen Tau Phenotypes") +
  ylab("Count") +
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_continuous(breaks = c(round(seq(-.6, .2, .1), 1)))

plt.nu <- ggplot(ps_gene[,c("Gene", "Nu.Avg", "Screen")], aes(Nu.Avg)) +
  geom_histogram(fill = "#999999") +
  geom_histogram(data = ps_gene[ps_gene$Screen,c("Gene", "Nu.Avg", "Screen")], 
                 aes(Nu.Avg), fill = "darkgreen") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Primary Screen Nu Phenotypes") +
  ylab("Count") +
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(1, 10, 100, 1000, 10000))

ggsave("../FIGURES/FIGURE-1/figure_1d_primary_screen_gamma_phenotype_distribution.svg", plt.gamma,
       width = 5, height = 5, dpi = 300)
ggsave("../FIGURES/FIGURE-1/figure_1d_primary_screen_tau_phenotype_distribution.svg", plt.tau,
       width = 5, height = 5, dpi = 300)
ggsave("../FIGURES/primary_screen_nu_phenotype_distribution.svg", plt.nu,
       width = 5, height = 5, dpi = 300)
