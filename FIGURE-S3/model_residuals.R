library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)

gphenos <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")
sphenos <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")

gp2 <- inner_join(gphenos[,c("sgRNA.id", "query", "single", "Expected", "Gamma.OI.Avg")], 
                  sphenos[,c("sgRNA.ID", "Gamma.Avg")], 
                  by = c("query" = "sgRNA.ID"))
gp2$additive <- gp2$single + gp2$Gamma.Avg

gp2 <- gp2[,c(1,2,3,6,4,7,5)]
colnames(gp2) <- c("a", "b", "ag", "bg", "quad", "add", "obs")

gp2$quad.gi <- gp2$obs - gp2$quad
gp2$add.gi <- gp2$obs - gp2$add

gp3 <- gp2 %>% 
  group_by(a) %>%
  summarise(mean.add = mean(add.gi),
            mean.quad = mean(quad.gi))
gp3 <- inner_join(gp3, 
                  sphenos[,c("sgRNA.ID", "Gamma.Avg")], 
                  by = c("a" = "sgRNA.ID"))

gp3.melt <- melt(gp3, id.vars = c("a","Gamma.Avg"))

plt <- ggplot(gp3.melt[nrow(gp3.melt):1,], aes(Gamma.Avg, value)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(aes(col = variable), alpha = .7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  xlab("Gene Gamma Phenotype") + ylab("Average Residual") +
  annotate(geom = "text", x = 0.01, y = .0625, hjust = 0,
           label = paste("r =",round(cor(gp3$mean.add, gp3$Gamma.Avg), 3))) +
  annotate(geom = "text", x = 0.01, y = .0575, hjust = 0,
           label = paste("r =",round(cor(gp3$mean.quad, gp3$Gamma.Avg), 3))) +
  scale_color_manual(values = c("black", "#999999"),
                     name = "Model",
                     labels = c("Additive", "Quadratic"))

ggsave("../FIGURES/TechnicalNote/model_residuals.svg", plt,
       height = 6, width =6, device = "svg")

