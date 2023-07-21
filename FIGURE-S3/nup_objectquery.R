gamma.r1.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_r1.txt")
gamma.r2.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_r2.txt")
gamma.avg.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")

phenos <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
id.map <- read_tsv("../DATA/genecombination_id_map.txt")

all.gamma <- inner_join(gamma.r1.all[,c(1:10,13,16)], 
                        gamma.r2.all[,c(1:10,13,16)], 
                        by = colnames(gamma.r1.all)[1:10], 
                        suffix = c(".R1", ".R2"))
all.gamma <- inner_join(all.gamma, gamma.avg.all[,c(1:10,13,16)], 
                        by = colnames(gamma.r1.all)[1:10])

nup <- "NUP62_+_50432687.23-P1P2"

nup.cons <- all.gamma[all.gamma$sgRNA.id == nup | all.gamma$query == nup,]
nup.cons$nonnup <- nup.cons$query
nup.cons$nonnup[nup.cons$nonnup == nup] <- nup.cons$sgRNA.id[nup.cons$nonnup == nup]
nup.cons$isquery <- nup.cons$query == nup
nup.cons$Avg <- rowMeans(nup.cons[,c("Gamma.OI.R1", "Gamma.OI.R2")])


nup.orient <- dcast(nup.cons[,c("isquery", "nonnup", "GI.z")], nonnup ~ isquery)
nup.orient <- nup.orient[!is.na(nup.orient$`TRUE`) & !is.na(nup.orient$`FALSE`),]
nup.orient2 <- inner_join(nup.orient, phenos[,c("sgRNA.ID", "Gamma.Avg")], by = c("nonnup" = "sgRNA.ID"))
nup.orient2$Diff <- nup.orient2$Gamma.Avg - phenos$Gamma.Avg[phenos$sgRNA.ID == "NUP62_+_50432687.23-P1P2"]

nup.orient.plt <- ggplot(nup.orient2, aes(`TRUE`, `FALSE`)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .6) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 2.5/1) +
  xlab("Score from NUP62 Model") +
  ylab("Score from Other Model") +
  scale_y_continuous(breaks = seq(-8,10,2)) 

ggsave("../FIGURES/TechnicalNote/nup_objectquery.svg", nup.orient.plt,
       device = "svg", height = 6, width = 3)
