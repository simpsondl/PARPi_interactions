library(readr)
library(dplyr)
library(ggplot2)

make_gamma_validation_plot <- function(gene1, gene2, single, interactions, pheno, col1, col2, col3, dc1, dc2, dc3){
  df <- data.frame(gene = c(rep(gene1, 2), rep(gene2, 2)),
                   pheno = c(rep(mean(as.numeric(unlist(single[grepl(gene1, single$sgRNA.ID), pheno]))), 2),
                             rep(mean(as.numeric(unlist(single[grepl(gene2, single$sgRNA.ID), pheno]))), 2)),
                   pts = c(as.numeric(unlist(single[grepl(gene1, single$sgRNA.ID), pheno])),
                           as.numeric(unlist(single[grepl(gene2, single$sgRNA.ID), pheno]))))
  df[5,] <- c("Expected", 
              mean(c(interactions$Expected[grepl(gene1, interactions$sgRNA.id) & grepl(gene2, interactions$query)],
                     interactions$Expected[grepl(gene2, interactions$sgRNA.id) & grepl(gene1, interactions$query)])),
              NA)
  df <- rbind(df,
              data.frame(gene = "Observed",
                         pheno = mean(c(as.numeric(unlist(interactions[grepl(gene1, interactions$sgRNA.id) & grepl(gene2, interactions$query), pheno])),
                                        as.numeric(unlist(interactions[grepl(gene2, interactions$sgRNA.id) & grepl(gene1, interactions$query), pheno])))),
                         pts = c(as.numeric(unlist(interactions[grepl(gene1, interactions$sgRNA.id) & grepl(gene2, interactions$query), pheno])),
                                 as.numeric(unlist(interactions[grepl(gene2, interactions$sgRNA.id) & grepl(gene1, interactions$query), pheno])))))
  
  df <- unique(df)
  
  df$fillgene <- df$gene
  df$fillgene[df$fillgene == gene1] <- "A"
  df$fillgene[df$fillgene == gene2] <- "B"
  df$fillpt <- paste0(df$fillgene, 2)
  df$fillgene <- factor(df$fillgene, levels = c("A", "B", "Expected", "Observed"))
  df$fillpt <- factor(df$fillpt, levels = c( "A2", "B2", "Expected2", "Observed2"))
  df$gene <- factor(df$gene, levels = c(gene1, gene2, "Expected", "Observed"))
  df$pheno <- as.numeric(df$pheno)
  df$pts <- as.numeric(df$pts)
  
  plt <- ggplot(df, aes(x = gene)) + 
    geom_hline(yintercept = -.3, linetype = "dotted") +
    geom_hline(yintercept = -.2, linetype = "dotted") +
    geom_hline(yintercept = -.1, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_col(data = unique(df[,c("gene", "fillgene", "pheno")]), 
             aes(y = pheno, fill = fillgene)) +
    geom_point(data = df[!is.na(df$pts),],
               aes(y = pts, fill = fillpt),
               shape = 21, color = "black", size = 2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.position = "none",
          aspect.ratio = 2/1) +
    scale_fill_manual(values = c(col1, dc1, col2, dc2, "#cccccc", col3, dc3 )) +
    xlab("") + ylab("Interaction Screen Gamma Phenotype") +
    scale_y_continuous(breaks = seq(-.3, 0, .1), limits = c(-.32, 0))
}

phenos <- read_tsv("../DATA/raw_phenotypes.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
ints <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")

phenos.sgc <- phenos %>% group_by(GuideCombinationID) %>%
  mutate(Gamma = mean(Gamma.Avg)) %>% select(FirstPosition:Orientation, 
                                             GuideCombinationID:PseudogeneCombinationID,
                                             Gamma) %>%
  filter(Orientation == "AB") %>%
  distinct()

single.pheno$Gamma.OI.Avg <- single.pheno$Gamma.Avg


# make_validation_plot <- function(guide1, guide2, g1label, g2label, single, interactions, pheno, col1, col2){
#   df <- data.frame(gene = c(g1label, g2label),
#                    pheno = c(as.numeric(single[single$sgRNA.ID == guide1 ,pheno]),
#                              as.numeric(single[single$sgRNA.ID == guide2 ,pheno])))
#   df[3,] <- c("Expected", mean(c(interactions$Expected[interactions$sgRNA.id == guide1 & interactions$query == guide2],
#                                interactions$Expected[interactions$sgRNA.id == guide2 & interactions$query == guide1])))
#   df[4,] <- c("Observed", mean(c(as.numeric(interactions[interactions$sgRNA.id == guide1 & interactions$query == guide2, pheno]),
#                                as.numeric(interactions[interactions$sgRNA.id == guide2 & interactions$query == guide1, pheno]))))
#   
#   df$gene <- factor(df$gene, levels = c(g1label, g2label, "Expected", "Observed"))
#   df$pheno <- as.numeric(df$pheno)
#   
#   plt <- ggplot(df, aes(gene, pheno, fill = gene)) + 
#     geom_col() + 
#     theme_bw() +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.grid = element_blank(),
#           axis.ticks.length = unit(.2, "cm"),
#           panel.border = element_blank(),
#           axis.line = element_line(),
#           legend.position = "none") +
#     scale_fill_manual(values = c(col1, col2, "#7e7f81", "black")) +
#     xlab("") + ylab("") 
# }
# 
# make_gene_validation_plot <- function(gene1, gene2, g1label, g2label, single, interactions, pheno, col1, col2){
#   df <- data.frame(gene = c(g1label, g2label),
#                    pheno = c(mean(as.numeric(unlist(single[grepl(gene1, single$sgRNA.ID), pheno]))),
#                              mean(as.numeric(unlist(single[grepl(gene2, single$sgRNA.ID), pheno])))))
#   df[3,] <- c("Expected", mean(c(interactions$Expected[grepl(gene1, interactions$sgRNA.id) & grepl(gene2, interactions$query)],
#                                  interactions$Expected[grepl(gene2, interactions$sgRNA.id) & grepl(gene1, interactions$query)])))
#   df[4,] <- c("Observed",mean(c(as.numeric(unlist(interactions[grepl(gene1, interactions$sgRNA.id) & grepl(gene2, interactions$query), pheno])),
#                                 as.numeric(unlist(interactions[grepl(gene2, interactions$sgRNA.id) & grepl(gene1, interactions$query), pheno])))))
#   
#   df$gene <- factor(df$gene, levels = c(g1label, g2label, "Expected", "Observed"))
#   df$pheno <- as.numeric(df$pheno)
#   
#   plt <- ggplot(df, aes(gene, pheno, fill = gene)) + 
#     geom_col() + 
#     theme_bw() +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.grid = element_blank(),
#           axis.ticks.length = unit(.2, "cm"),
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           panel.border = element_blank(),
#           axis.line = element_line(),
#           legend.position = "none",
#           aspect.ratio = 2/1) +
#     scale_fill_manual(values = c(col1, col2, "#7e7f81", "black")) +
#     xlab("") + ylab("Interaction Screen Gamma Phenotype") +
#     scale_y_continuous(breaks = seq(-.3, 0, .05), limits = c(-.3, 0))
# }


plt2 <- make_gamma_validation_plot("MND1", "RAD54L",
                                   single.pheno, ints,
                                   "Gamma.OI.Avg", 
                                   "#ddcc77", "#117733", "black",
                                   "#b59322", "#519751", "#666666")

plt3 <- make_gamma_validation_plot("PSMC3IP", "RAD54L",
                                   single.pheno, ints,
                                   "Gamma.OI.Avg", 
                                   "#CC6677", "#117733", "black",
                                   "#ee99aa", "#519751", "#666666")

plt4 <- make_gamma_validation_plot("MND1", "PSMC3IP",
                                   single.pheno, ints,
                                   "Gamma.OI.Avg", 
                                   "#ddcc77", "#CC6677", "black",
                                   "#b59322", "#ee99aa", "#666666")

ggsave("../FIGURES/FIGURE-3/mnd1_rad54l_gammascreen_phenotypes.svg", plt2,
       device = "svg", height = 4, width = 2)
ggsave("../FIGURES/FIGURE-3/psmc3ip_rad54l_gammascreen_phenotypes.svg", plt3,
       device = "svg", height = 4, width = 2)
ggsave("../FIGURES/FIGURE-3/mnd1_psmc3ip_gammascreen_phenotypes.svg", plt4,
       device = "svg", height = 4, width = 2)
