library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load functions
source("helper_functions.R")

# Load data
annotations <- read_tsv("../DATA/External_Data/screen_annotations_w_gene_inclusion_info.txt")
# Load primary screen data and calculate phenotypes
primary_counts <- read_delim("../DATA/raw_counts_primary_screen.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE, skip = 4)
colnames(primary_counts) <- c("sgRNA.ID","NIRAP.R1","NIRAP.R2","T0.R1","T0.R2", "DMSO.R1", "DMSO.R2")

mageck_pvalues <- read_tsv("../DATA/External_Data/mageck_gene_results_primary_screen.txt")
mageck_pvalues$Gene2 <- gsub("-.*", "", mageck_pvalues$Gene)
mageck_pvalues$twosided[mageck_pvalues$twosided > 1] <- apply(mageck_pvalues[mageck_pvalues$twosided > 1,c("neg|p-value", "pos|p-value")], 
                                                              1, FUN = function(x) 2*(1 - min(x)))
mageck <- mageck_pvalues %>% group_by(Gene2) %>% summarise(twosided = min(twosided))

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

ps_gene_phenos <- inner_join(ps_gene_phenos, mageck, by = c("Gene" = "Gene2"))

annotations$Highlight <- NA
annotations$Highlight[annotations$PrimaryScreenHit & annotations$N.LiteratureHit == 0] <- "Interaction library - unique to primary screen"
annotations$Highlight[annotations$N.LiteratureHit > 0] <- "Interaction library - any literature screen"
annotations$Highlight[is.na(annotations$Highlight)] <- "Interaction library - included"

ps_gene_phenos <- left_join(ps_gene_phenos, annotations[,c("MapGeneName", "Highlight")], by = c("Gene" = "MapGeneName"))
ps_gene_phenos$Highlight[is.na(ps_gene_phenos$Highlight)] <- "Not present"

ps_gene_phenos$Highlight <- factor(ps_gene_phenos$Highlight, levels = c("Not present",
                                                                        "Interaction library - included",
                                                                        "Interaction library - any literature screen",
                                                                        "Interaction library - unique to primary screen"))

ps_gene_phenos <- ps_gene_phenos[order(ps_gene_phenos$Highlight),]

ps_gene_phenos$Lbl <- NA
ps_gene_phenos$Lbl[(ps_gene_phenos$Rho.Avg < -1.02 | ps_gene_phenos$Rho.Avg > .4) &
                     ps_gene_phenos$twosided == min(ps_gene_phenos$twosided)] <- ps_gene_phenos$Gene[(ps_gene_phenos$Rho.Avg < -1.02 | ps_gene_phenos$Rho.Avg > .4) &
                                                                                                       ps_gene_phenos$twosided == min(ps_gene_phenos$twosided)]
#ps_gene_phenos$Lbl[ps_gene_phenos$Highlight == "Interaction library - unique to primary screen" & 
#                     ps_gene_phenos$twosided == min(ps_gene_phenos$twosided)] <- ps_gene_phenos$Gene[ps_gene_phenos$Highlight == "Interaction library - unique to primary screen" & ps_gene_phenos$twosided == min(ps_gene_phenos$twosided)]


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

plt <- ggplot(ps_gene_phenos, aes(Rho.Avg, -log10(twosided))) + 
  geom_point(data = ps_gene_phenos[ps_gene_phenos$Highlight == "Not present",], 
             alpha = 1, size = .7, col = "#bbbbbb") +
  geom_point(data = ps_gene_phenos[ps_gene_phenos$Highlight != "Not present",], 
             aes(col = Highlight), size = .9, alpha = .6) +
  theme_bw() +
  geom_text_repel(data = ps_gene_phenos[!is.na(ps_gene_phenos$Lbl),], 
                  aes(label = Lbl), min.segment.length = 0) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm")) +
  scale_color_manual(values = c("black", "dodgerblue", "#BB5566")) +
  xlim(c(-1.6,1.6)) +
  ylim(c(-.01, 7)) +
  xlab("Rho Phenotype") + labs(y = expression("-log"["10"] ~ "(p-value)")) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

ggsave("../FIGURES/FIGURE-1/primary_screen_volcano_plot.svg", plt, device = "svg", height = 6, width = 8, dpi = 300)
  