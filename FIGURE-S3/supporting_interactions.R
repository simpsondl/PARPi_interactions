library(readr)
library(ggplot2)
library(dplyr)

gamma.r1.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_r1.txt")
gamma.r2.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_r2.txt")
gamma.avg.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")

phenos <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
id.map <- read_tsv("../DATA/genecombination_id_map.txt")

all.gamma <- inner_join(gamma.r1.all[,c(1:10,13,16)], 
                        gamma.r2.all[,c(1:10,13,16)], 
                        by = colnames(gamma.r1.all)[1:10], 
                        suffix = c(".R1", ".R2"))
all.gamma <- inner_join(all.gamma, gamma.avg.all[,c(1:10,13,16)], 
                        by = colnames(gamma.r1.all)[1:10])

gamma.gcs <- all.gamma %>% group_by(GeneCombinationID) %>% 
  summarise(GI.R1 = mean(GI.z.R1), GI.R2 = mean(GI.z.R2), GI.Avg = mean(GI.z), N = n())
gamma.gcs <- inner_join(gamma.gcs, id.map)

gc.var <- gamma.gcs %>% 
  group_by(N) %>% 
  filter(N <= 8) %>% 
  summarise(mean = mean(GI.Avg), var = var(GI.Avg), sig = sum(abs(GI.Avg) >= 1.5274), Count = n())
gc.var$frac <- round(gc.var$sig/gc.var$Count * 100)
gamma.gcs2 <- inner_join(gamma.gcs, gc.var)
gc.counts <- as.data.frame(table(gamma.gcs2$N))

supp.plt <- ggplot(gamma.gcs2, aes(as.factor(N), GI.Avg)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_hline(yintercept = -1.5274, alpha = .5, linetype = "dashed") +
  geom_hline(yintercept = 1.5274, alpha = .5, linetype = "dashed") +
  geom_violin(aes(fill = var)) +
  geom_text(data = gc.var, aes(x = as.factor(N), y = -16, label = Count), size = 3) +
  geom_text(data = gc.var, aes(x = as.factor(N), y = -17, label = paste0("(", frac, "%)")), size = 2.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.background = element_blank()) +
  xlab("Number of Supporting Interactions") +
  ylab("Interaction Score") + 
  scale_fill_viridis_c(name = "Variance")

ggsave("../FIGURES/TechnicalNote/all_supporting.pdf", supp.plt, device = "pdf",
       height = 6, width = 6, dpi = 300)
ggsave("../FIGURES/TechnicalNote/all_supporting.svg", supp.plt, device = "svg",
       height = 6, width = 6)
