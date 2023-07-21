library(readr)
library(dplyr)
library(ggplot2)

phenos <- read_tsv("../DATA/raw_rho_phenotypes.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_rho_phenotypes.txt")

clusts <- read_tsv("../DATA/heatmap_clusters_power4power6.txt")
tauclusts <- read_tsv("../DATA/tau_clusters_p4.txt")

# DEFINED BY LITERATURE REVIEW
fancs <- c("FANCB", "FANCA", "FANCG", "FANCF", "FANCM", "C19orf40",  
           "FANCE", "FANCC", "FANCD2", "FANCI", "UBE2T",
           "BRIP1", "RAD18")
brca1a <- c("UIMC1", "BRCC3", "BABAM1", "FAM175A", "BRCA1")
buffs <- c("TP53BP1", "RIF1", "MAD2L2", "FAM35A", "C20orf196", "CCAR1",  "CTC1", 
           "OBFC1", "TEN1", "ATM")
buffs2 <- c(buffs, "PAXIP1", "POLA1", "POLA2", "PRIM1", "PRIM2")

ssi <- c(unique(c(fancs, brca1a, buffs)), "AUNIP")

single.pheno$gene <- gsub("_.*", "", single.pheno$sgRNA.ID)
single.gene <- single.pheno %>% group_by(gene) %>%
  mutate(Rho = mean(Rho.Avg)) %>%
  filter(gene %in% ssi)

single.gene$group <- NA
single.gene$group[single.gene$gene %in% setdiff(fancs, "BRCA1")] <- "Fanconi"
single.gene$group[is.na(single.gene$group)] <- ifelse(single.gene$gene[is.na(single.gene$group)] %in% brca1a, "BRCA1A", "CST+")
single.gene$group[single.gene$gene == "AUNIP"] <- "AUNIP"

single.gene <- single.gene[order(single.gene$Rho),]

single.gene$group <- factor(single.gene$group, levels= c("AUNIP", "BRCA1A", "CST+", "Fanconi"))
single.gene$gene <- factor(single.gene$gene, levels = c("AUNIP", 
                                                        "BRCA1", "FAM175A", "UIMC1", "BRCC3", "BABAM1",
                                                        "RIF1", "ATM", "MAD2L2", "TEN1",
                                                        "OBFC1", "FAM35A", "C20orf196", "CTC1", 
                                                        "CCAR1", "TP53BP1",
                                                        'FANCB', 'FANCA', 'C19orf40', 'FANCD2', 'FANCI', 'UBE2T', 
                                                        'FANCM', 'FANCF', 'FANCG', 'FANCE', 'RAD18', 'FANCC', 'BRIP1'))

single.gene$fill <- paste0(single.gene$group, "2")

plt <- ggplot(single.gene, aes(gene)) +
  facet_grid(~ group, scales = "free_x", space = "free") +
  geom_hline(yintercept = -.9, linetype = "dotted") +
  geom_hline(yintercept = -.75, linetype = "dotted") +
  geom_hline(yintercept = -.6, linetype = "dotted") +
  geom_hline(yintercept = -.45, linetype = "dotted") +
  geom_hline(yintercept = -.3, linetype = "dotted") +
  geom_hline(yintercept = -.15, linetype = "dotted") +
  geom_hline(yintercept = .15, linetype = "dotted") +
  geom_col(data = unique(single.gene[,c("gene", "Rho", "group")]), 
           aes(y= Rho, fill = group)) +
  geom_point(aes(y = Rho.Avg, fill = fill), 
             col = "black", shape = 21) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_manual(values = c("black", "#666666", "#117733", "#519751",
                               "#DDCC77", "#b59322",
                               "#CC6677", "#ee99aa")) +
  scale_y_continuous(breaks = seq(-.9, 0, .3),
                     labels = c(-.9, -.6, -.3, 0)) +
  xlab("") + ylab("Screen Rho Phenotype")

ggsave("../FIGURES/FIGuRE-6/aunip_resection_phenotype_barplot_072023.svg", plt,
       device = "svg", height = 3, width = 6)
