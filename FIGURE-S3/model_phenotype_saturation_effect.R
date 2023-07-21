library(dplyr)
library(readr)

gquad <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_gamma_oi_avg.txt")
tquad <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_tau_oi_avg.txt")
glin <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_linear_gamma_oi_avg.txt")
tlin <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_linear_tau_oi_avg.txt")

gquadst <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_statistics_gamma_oi_avg.txt")
tquadst <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_statistics_tau_oi_avg.txt")
glinst <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_statistics_linear_gamma_oi_avg.txt")
tlinst <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_statistics_linear_tau_oi_avg.txt")

pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
horlbeck <- read_tsv("F:/Downloads/NIHMS973971-supplement-11.txt")

horlbeck$guide1 <- gsub("\\+\\+.*","", horlbeck$...1)
horlbeck$guide2 <- gsub(".*\\+\\+","", horlbeck$...1)
horlbeck$gene1 <- gsub("_.*", "", horlbeck$guide1)
horlbeck$gene2 <- gsub("_.*", "", horlbeck$guide2)
horlbeck$Category <- "X+Y"
horlbeck$Category[horlbeck$gene1 == horlbeck$gene2] <- "X+X"
horlbeck$Category[horlbeck$gene1 == "negative" |  horlbeck$gene2 == "negative"] <- "X+NT"
horlbeck$Category[horlbeck$gene1 == "negative" &  horlbeck$gene2 == "negative"] <- "NT+NT"

horlbeck2 <- horlbeck[3:nrow(horlbeck),c(1,5:ncol(horlbeck))]
colnames(horlbeck2)[1:4] <- c("GuideCombinationName", "Gamma.R1", "Gamma.R2", "Gamma.Avg")

horlbeck2$Gamma.R1 <- as.numeric(horlbeck2$Gamma.R1)
horlbeck2$Gamma.R2 <- as.numeric(horlbeck2$Gamma.R2)
horlbeck2$Gamma.Avg <- as.numeric(horlbeck2$Gamma.Avg)


h.pheno <- data.frame("sgRNA.ID" = unique(horlbeck2$guide1), 
                      "Gamma.R1" = 0, 
                      "Gamma.R2" = 0,
                      "Gamma.Avg" = 0)

for(i in 1:nrow(h.pheno)){
  # Handle non-targeting case
  if(grepl("negative",h.pheno$sgRNA.ID[i])){
    # Extract all non-targeting guide combinations with desired non-targeting guide in position A or B
    # If this is not handled explicitly, will aggregate all NT guides into a glob
    tmp <- horlbeck2[(horlbeck2$guide1 == h.pheno$sgRNA.ID[i] &
                        horlbeck2$guide2 %in% unique(horlbeck2$guide1[horlbeck2$Category == "NT+NT"])) | 
                       (horlbeck2$guide2 == h.pheno$sgRNA.ID[i] &
                          horlbeck2$guide1 %in% unique(horlbeck2$guide1[horlbeck2$Category == "NT+NT"])),]
    h.pheno$Gamma.R1[i] <- mean(tmp$Gamma.R1)
    h.pheno$Gamma.R2[i] <- mean(tmp$Gamma.R2)
    h.pheno$Gamma.Avg[i] <- mean(tmp$Gamma.Avg)
    h.pheno$N[i] <- nrow(tmp)
    
  }
  # Handle targeting case
  else {
    # Extract all combinations with desired targeting guide in position A or B and non-targeting guide in other
    # Easier to do this since we can just grab all non-targeting at once
    tmp <- horlbeck2[(horlbeck2$guide1 == h.pheno$sgRNA.ID[i] | horlbeck2$guide2 == h.pheno$sgRNA.ID[i]) & 
                       (horlbeck2$guide1 %in% unique(horlbeck2$guide1[horlbeck2$Category == "NT+NT"]) | 
                          horlbeck2$guide2 %in% unique(horlbeck2$guide1[horlbeck2$Category == "NT+NT"])),]
    h.pheno$Gamma.R1[i] <- mean(tmp$Gamma.R1)
    h.pheno$Gamma.R2[i] <- mean(tmp$Gamma.R2)
    h.pheno$Gamma.Avg[i] <- mean(tmp$Gamma.Avg)
    h.pheno$N[i] <- nrow(tmp)
    
  }
}

h.pheno2 <- h.pheno[!is.na(h.pheno$Gamma.Avg),]

gcomp <- inner_join(glinst, 
                    gquadst,
                    by = "sgRNA.ID", suffix = c(".lin", ".quad"))
gcomp$F <- (gcomp$r.squared.quad - gcomp$r.squared.lin) * gcomp$df.residual.quad/(1 - gcomp$r.squared.quad)
gcomp$pf <- pf(gcomp$F, gcomp$df.residual.lin, gcomp$df.residual.quad, lower.tail = FALSE)
glinguides <- gcomp$sgRNA.ID[gcomp$pf >= .05/nrow(gcomp)]

glin2 <- inner_join(glin, pheno[,c("sgRNA.ID", "Gamma.Avg")]) %>%
  filter(sgRNA.ID %in% glinguides)

ggplot(glin2[glin2$term == "single",], 
       aes(Gamma.Avg, estimate)) +
  geom_point() +
  geom_smooth(method = "lm")

hquad <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_horlbeckquad_gamma_oi_avg.txt")
hlin <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_horlbecklinear_gamma_oi_avg.txt")

hquadst <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_statistics_horlbeckquad_gamma_oi_avg.txt")
hlinst <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_statistics_horlbecklinear_gamma_oi_avg.txt")

hcomp <- inner_join(hlinst[,c("sgRNA.ID", "r.squared", "df.residual")],
                    hquadst[,c("sgRNA.ID", "r.squared", "df.residual")],
                    by = "sgRNA.ID", suffix = c(".lin", ".quad"))
hcomp$F <- (hcomp$r.squared.quad - hcomp$r.squared.lin) * hcomp$df.residual.quad/(1 - hcomp$r.squared.quad)
hcomp$pf <- pf(hcomp$F, hcomp$df.residual.lin, hcomp$df.residual.quad, lower.tail = FALSE)
hlinguides <- hcomp$sgRNA.ID[hcomp$pf >= .05/nrow(hcomp)]

hlin2 <- inner_join(hlin, h.pheno2[,c("sgRNA.ID", "Gamma.Avg")]) %>%
  filter(sgRNA.ID %in% hlinguides)

ggplot(hlin2[hlin2$term == "single",], 
       aes(Gamma.Avg, estimate)) +
  geom_point() +
  geom_smooth(method = "lm")

hfit <- lm(estimate ~ Gamma.Avg, data = hlin2)

glin2$Screen <- "Simpson"
hlin2$Screen <- "Horlbeck"

cmp <- rbind(glin2[glin2$term == "single", c("sgRNA.ID", "estimate", "Gamma.Avg", "Screen")],
             hlin2[hlin2$term == "single", c("sgRNA.ID", "estimate", "Gamma.Avg", "Screen")])

plt <- ggplot(cmp, aes(Gamma.Avg, estimate, col = Screen)) + 
  geom_hline(yintercept = 1, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  xlab("Guide Gamma Phenotype") + ylab("Slope of Linear Fit") +
  scale_color_manual(values = c("#bbbbbb", "black")) +
  annotate(geom = "text", x = -.7, y = 1.15,
           label = paste0("rs = ", round(cor(glin2$estimate[glin2$term == "single"],
                                             glin2$Gamma.Avg[glin2$term == "single"]), 3))) +
  annotate(geom = "text", x = -.7, y = 1.05,
           label = paste0("rh = ", round(cor(hlin2$estimate[hlin2$term == "single"],
                                             hlin2$Gamma.Avg[hlin2$term == "single"]), 3)))

ggsave("../FIGURES/TechnicalNote/modification_of_additive.svg", plt,
       height = 6, width = 6, device = "svg")
