library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

model <- read_tsv("../DATA/Interaction_Scores/Gamma.OI.Avg/BRCA1_-_41277328.23-P1P2.txt")
model_coeffs <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_gamma_oi_avg.txt")
model_coeffs <- model_coeffs[model_coeffs$sgRNA.ID == "BRCA1_-_41277328.23-P1P2",]

model$GI.capped <- model$GI.z
model$GI.capped[model$GI.capped >= 3.5] <- 3.5
model$GI.capped[model$GI.capped <= -3.5] <- -3.5

model$lb <- NA
model$lb[model$sgRNA.id == "TP53BP1_+_43785309.23-P1"] <- "TP53BP1_1"
model$lb[model$sgRNA.id == "TP53BP1_+_43785094.23-P1"] <- "TP53BP1_2"
model$lb[model$sgRNA.id == "BARD1_+_215674399.23-P1P2"] <- "BARD1_1"
model$lb[model$sgRNA.id == "BARD1_+_215674388.23-P1P2"] <- "BARD1_2"
model$lb[model$sgRNA.id == "C19orf40_-_33463200.23-P1P2"] <- "C19orf40_1"
model$lb[model$sgRNA.id == "C19orf40_-_33463205.23-P1P2"] <- "C19orf40_2"
model$lb[model$sgRNA.id == "FANCA_-_89883018.23-P1P2"] <- "FANCA_1"
model$lb[model$sgRNA.id == "FANCA_+_89883007.23-P1P2"] <- "FANCA_2"
model$lb[model$sgRNA.id == "FANCB_-_14891130.23-P1P2"] <- "FANCB_1"
model$lb[model$sgRNA.id == "FANCB_-_14891133.23-P1P2"] <- "FANCB_2"
model$lb[model$sgRNA.id == "FANCD2_-_10068143.23-P1P2"] <- "FANCD2_1"
model$lb[model$sgRNA.id == "FANCD2_-_10068176.23-P1P2"] <- "FANCD2_2"
model$lb[model$sgRNA.id == "FANCF_-_22647359.23-P1P2"] <- "FANCF_1"
model$lb[model$sgRNA.id == "FANCF_-_22646985.23-P1P2"] <- "FANCF_2"
model$lb[model$sgRNA.id == "FANCG_+_35079731.23-P1P2"] <- "FANCG_1"
model$lb[model$sgRNA.id == "FANCG_-_35079943.23-P1P2"] <- "FANCG_2"
model$lb[model$sgRNA.id == "FANCM_+_45605620.23-P1P2"] <- "FANCM_1"
model$lb[model$sgRNA.id == "FANCM_-_45605167.23-P1P2"] <- "FANCM_2"
model$lb[model$sgRNA.id == "UBE2T_-_202311067.23-P1P2"] <- "UBE2T_1"
model$lb[model$sgRNA.id == "UBE2T_-_202310883.23-P1P2"] <- "UBE2T_2"



model_line <- data.frame(single = seq(-.37, .07, .01), pred = NA)
model_line$pred <- model_coeffs$estimate[1] + model_coeffs$estimate[2] * model_line$single + model_coeffs$estimate[3] * model_line$single^2
model_line$lb <- NA

plt <- ggplot(model, aes(single, Gamma.OI.Avg, col = GI.capped, label = lb)) + 
        geom_line(data = model_line, aes(single, pred), col = "red", size = 2) +
        geom_point() + 
        geom_point(data = model[!is.na(model$lb),], shape = 1, size = 3, col = "black") +
        geom_text_repel(size = 3, max.overlaps = Inf, min.segment.length = 0) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.grid = element_blank(),
              axis.ticks.length = unit(.2, "cm")) +
        xlab("Single sgRNA Gamma") +
        ylab("Combination Gamma") +
        scale_color_gradient2(name = "IS",
                              low = "dodgerblue", 
                              mid = "black", 
                              high = "darkorange1", 
                              midpoint = mean(model$Gamma.OI.Avg[model$Control]),
                              breaks = c(-3, -1.5, 0, 1.5, 3))

plt.nolab <- ggplot(model, aes(single, Gamma.OI.Avg, col = GI.capped)) + 
  geom_line(data = model_line, aes(single, pred), col = "red", size = 2) +
  geom_point() + 
  geom_point(data = model[!is.na(model$lb),], shape = 1, size = 3, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  xlab("") +
  ylab("") +
  scale_color_gradient2(name = "IS",
                        low = "dodgerblue", 
                        mid = "black", 
                        high = "darkorange1", 
                        midpoint = mean(model$Gamma.OI.Avg[model$Control]),
                        breaks = c(-3, -1.5, 0, 1.5, 3)) 

ggsave("../FIGURES/FIGURE-2/brca1_model_scatterplot_labels_withlines_circled.svg", plt, device = "svg", height = 6, width = 7, dpi = 300)
ggsave("../FIGURES/FIGURE-2/brca1_model_scatterplot_nolabels.svg", plt.nolab, device = "svg", height = 6, width = 6, dpi = 300)
