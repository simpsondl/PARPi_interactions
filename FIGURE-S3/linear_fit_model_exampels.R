library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)

phc <- read_tsv("../DATA/Interaction_Scores/Gamma.OI.Avg/PHC2_-_33815107.23-P1P2.txt")
znf <- read_tsv("../DATA/Interaction_Scores/Gamma.OI.Avg/ZNF574_-_42580338.23-P1.txt")
emc <- read_tsv("../DATA/Interaction_Scores/Gamma.OI.Avg/EMC2_+_109456085.23-P1P2.txt")
model_coeffs <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_linear_gamma_oi_avg.txt")
model_coeffs <- model_coeffs[model_coeffs$sgRNA.ID %in% c("PHC2_-_33815107.23-P1P2",
                                                          "ZNF574_-_42580338.23-P1",
                                                          "EMC2_+_109456085.23-P1P2"),]
phc$gene <- 'PHC2'
znf$gene <- 'ZNF574'
emc$gene <- "EMC2"

models <- rbind(phc[,c("gene", "single", "Gamma.OI.Avg")],
                znf[,c("gene", "single", "Gamma.OI.Avg")],
                emc[,c("gene", "single", "Gamma.OI.Avg")])

models$gene <- factor(models$gene, levels = c("PHC2", "EMC2", "ZNF574"))

model_line <- data.frame(single = seq(-.37, .07, .01), 
                         PHC2 = NA,
                         ZNF574 = NA,
                         EMC2 = NA)

model_line$PHC2 <- model_coeffs$estimate[3] + model_coeffs$estimate[4] * model_line$single
model_line$ZNF574 <- model_coeffs$estimate[5] + model_coeffs$estimate[6] * model_line$single
model_line$EMC2 <- model_coeffs$estimate[1] + model_coeffs$estimate[2] * model_line$single

mline_melt <- melt(model_line, id.vars = "single")
mline_melt$variable <- factor(mline_melt$variable, levels = c("PHC2", "EMC2", "ZNF574"))

plt <- ggplot(models, aes(single, Gamma.OI.Avg, col = gene)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_line(data = mline_melt, aes(single, value, col = variable), size = 1.2) +
  geom_point(size = .5, alpha = .7) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Single sgRNA Gamma") +
  ylab("Combinatorial Gamma") +
  scale_color_manual(values = c("#DDAA33", "#BB5566", "#004488"),
                     name = "Gene")

ggsave("../FIGURES/TechnicalNote/linear_fit_model_examples.svg", plt,
       device = "svg", height = 5, width = 5)
