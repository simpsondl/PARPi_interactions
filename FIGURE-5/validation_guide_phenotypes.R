library(readr)
library(dplyr)
library(ggplot2)

phenos <- read_tsv("../DATA/raw_phenotypes.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")
ints <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")

phenos.sgc <- phenos %>% group_by(GuideCombinationID) %>%
  mutate(Gamma = mean(Gamma.Avg)) %>% select(FirstPosition:Orientation, 
                                             GuideCombinationID:PseudogeneCombinationID,
                                             Gamma) %>%
  filter(Orientation == "AB") %>%
  distinct()

parp1.2 <- "PARP1_+_226595771.23-P1P2"
parp2.1 <- "PARP2_-_20811828.23-P1P2"
parp2.2 <- "PARP2_-_20811834.23-P1P2"
rnaseh2c.1 <- "RNASEH2C_-_65488303.23-P1P2"
rnaseh2c.2 <- "RNASEH2C_+_65488196.23-P1P2"

df <- data.frame(gc = c("PARP1 (2)", "RNASEH2C (1)"), gamma = c(single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == parp1.2],
                                                        single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == rnaseh2c.1]))
df[3,] <- c("Expected", mean(ints$Expected[ints$sgRNA.id == parp1.2 & ints$query == rnaseh2c.1],
                             ints$Expected[ints$sgRNA.id == rnaseh2c.1 & ints$query == parp1.2]))
df[4,] <- c("Observed", phenos.sgc$Gamma[phenos.sgc$FirstPosition == parp1.2 & phenos.sgc$SecondPosition == rnaseh2c.1])                 

df$gc <- factor(df$gc, levels = c("PARP1 (2)", "RNASEH2C (1)", "Expected", "Observed"))
df$gamma <- as.numeric(df$gamma)

parp12.rnaseh2c1 <- ggplot(df, aes(gc, gamma, fill = gc)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none") +
  scale_fill_manual(values = c("#ddcc77", "#CC6677", "#7e7f81", "black")) +
  xlab("")

ggsave("../FIGURES/FIGURE-5/parp1-2_rnaseh2c-1.svg", parp12.rnaseh2c1,
       device = "svg", height = 4, width = 4)

#####################################

df2 <- data.frame(gc = c("PARP2 (1)", "RNASEH2C (1)"), gamma = c(single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == parp2.1],
                                                                single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == rnaseh2c.1]))
df2[3,] <- c("Expected", mean(ints$Expected[ints$sgRNA.id == parp2.1 & ints$query == rnaseh2c.1],
                             ints$Expected[ints$sgRNA.id == rnaseh2c.1 & ints$query == parp2.1]))
df2[4,] <- c("Observed", phenos.sgc$Gamma[phenos.sgc$FirstPosition == parp2.1 & phenos.sgc$SecondPosition == rnaseh2c.1])                 

df2$gc <- factor(df2$gc, levels = c("PARP2 (1)", "RNASEH2C (1)", "Expected", "Observed"))
df2$gamma <- as.numeric(df2$gamma)

parp21.rnaseh2c1 <- ggplot(df2, aes(gc, gamma, fill = gc)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line()) +
  scale_fill_manual(values = c("#117733", "#CC6677", "#7e7f81", "black")) +
  xlab("")

ggsave("../FIGURES/FIGURE-5/parp2-1_rnaseh2c-1.svg", parp21.rnaseh2c1,
       device = "svg", height = 4, width = 4)

#####################################

df3 <- data.frame(gc = c("PARP1", "RNASEH2C (2)"), gamma = c(mean(single.pheno$Gamma.Avg[grepl("PARP1",single.pheno$sgRNA.ID)]),
                                                                 single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == rnaseh2c.2]))
df3[3,] <- c("Expected", mean(c(ints$Expected[grepl("PARP1",ints$sgRNA.id) & ints$query == rnaseh2c.2],
                              ints$Expected[ints$sgRNA.id == rnaseh2c.2 & grepl("PARP1",ints$query)])))
df3[4,] <- c("Observed", mean(phenos.sgc$Gamma[grepl("PARP1", phenos.sgc$FirstPosition) & phenos.sgc$SecondPosition == rnaseh2c.2]))                 

df3$gc <- factor(df3$gc, levels = c("PARP1", "RNASEH2C (2)", "Expected", "Observed"))
df3$gamma <- as.numeric(df3$gamma)

parp1g.rnaseh2c2 <- ggplot(df3, aes(gc, gamma, fill = gc)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("#45A1DA", "#F78C1E", "#7e7f81", "black")) +
  xlab("")

ggsave("../FIGURES/FIGURE-5/parp1-gene_rnaseh2c-2.svg", parp1g.rnaseh2c2,
       device = "svg", height = 4, width = 4)

#####################################

df4 <- data.frame(gc = c("PARP2 (2)", "RNASEH2C (2)"), gamma = c(single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == parp2.2],
                                                                 single.pheno$Gamma.Avg[single.pheno$sgRNA.ID == rnaseh2c.2]))
df4[3,] <- c("Expected", mean(ints$Expected[ints$sgRNA.id == parp2.2 & ints$query == rnaseh2c.2],
                              ints$Expected[ints$sgRNA.id == rnaseh2c.2 & ints$query == parp2.2]))
df4[4,] <- c("Observed", phenos.sgc$Gamma[phenos.sgc$FirstPosition == parp2.2 & phenos.sgc$SecondPosition == rnaseh2c.2])                 

df4$gc <- factor(df4$gc, levels = c("PARP2 (2)", "RNASEH2C (2)", "Expected", "Observed"))
df4$gamma <- as.numeric(df4$gamma)

parp22.rnaseh2c2 <- ggplot(df4, aes(gc, gamma, fill = gc)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("#a071b1", "#F78C1E", "#7e7f81", "black")) +
  xlab("")

ggsave("../FIGURES/FIGURE-5/parp2-2_rnaseh2c-2.svg", parp22.rnaseh2c2,
       device = "svg", height = 4, width = 4)

#####################################

df5 <- data.frame(gc = c("PARP1", "RNASEH2C"), gamma = c(mean(single.pheno$Gamma.Avg[grepl("PARP1",single.pheno$sgRNA.ID)]),
                                                             mean(single.pheno$Gamma.Avg[grepl("RNASEH2C",single.pheno$sgRNA.ID)])))
df5[3,] <- c("Expected", mean(c(ints$Expected[grepl("PARP1",ints$sgRNA.id) & grepl("RNASEH2C",ints$query)],
                                ints$Expected[grepl("RNASEH2C",ints$sgRNA.id) & grepl("PARP1",ints$query)])))
df5[4,] <- c("Observed", mean(phenos.sgc$Gamma[grepl("PARP1", phenos.sgc$FirstPosition) & grepl("RNASEH2C", phenos.sgc$SecondPosition)]))                 

df5$gc <- factor(df5$gc, levels = c("PARP1", "RNASEH2C", "Expected", "Observed"))
df5$gamma <- as.numeric(df5$gamma)

parp1g.rnaseh2cg <- ggplot(df5, aes(gc, gamma, fill = gc)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line()) +
  scale_fill_manual(values = c("#ddcc77", "#CC6677", "#7e7f81", "black")) +
  xlab("")

ggsave("../FIGURES/FIGURE-5/parp1-gene_rnaseh2c-gene.svg", parp1g.rnaseh2cg,
       device = "svg", height = 4, width = 4)


#####################################

df6 <- data.frame(gc = c("PARP2", "RNASEH2C"), gamma = c(mean(single.pheno$Gamma.Avg[grepl("PARP2",single.pheno$sgRNA.ID)]),
                                                         mean(single.pheno$Gamma.Avg[grepl("RNASEH2C",single.pheno$sgRNA.ID)])))
df6[3,] <- c("Expected", mean(c(ints$Expected[grepl("PARP2",ints$sgRNA.id) & grepl("RNASEH2C",ints$query)],
                                ints$Expected[grepl("RNASEH2C",ints$sgRNA.id) & grepl("PARP2",ints$query)])))
df6[4,] <- c("Observed", mean(phenos.sgc$Gamma[grepl("PARP2", phenos.sgc$FirstPosition) & grepl("RNASEH2C", phenos.sgc$SecondPosition)]))                 

df6$gc <- factor(df6$gc, levels = c("PARP2", "RNASEH2C", "Expected", "Observed"))
df6$gamma <- as.numeric(df6$gamma)

parp2g.rnaseh2cg <- ggplot(df6, aes(gc, gamma, fill = gc)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line()) +
  scale_fill_manual(values = c("#117733", "#CC6677", "#7e7f81", "black")) +
  xlab("")

ggsave("../FIGURES/FIGURE-5/parp2-gene_rnaseh2c-gene.svg", parp2g.rnaseh2cg,
       device = "svg", height = 4, width = 4)
