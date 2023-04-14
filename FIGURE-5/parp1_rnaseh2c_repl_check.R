library(readr)
library(dplyr)
library(ggplot2)

source("helper_functions.R")

parp1.1 <- "PARP1_+_226595744.23-P1P2"
parp1.2 <- "PARP1_+_226595771.23-P1P2"
rnasea.1 <- "RNASEH2A_-_12917521.23-P1P2"
rnasea.2 <- "RNASEH2A_+_12917430.23-P1P2"
rnaseb.1 <- "RNASEH2B_+_51483920.23-P1P2"
rnaseb.2 <- "RNASEH2B_+_51483924.23-P1P2"
rnasec.1 <- "RNASEH2C_-_65488303.23-P1P2"
rnasec.2 <- "RNASEH2C_+_65488196.23-P1P2"

# Import scores for BRCA2 models and all interaction scores
# Note that these files are subsets of the imports of the two files below
gamma.r1.parp1.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/PARP1_+_226595744.23-P1P2.txt")
gamma.r1.parp1.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/PARP1_+_226595771.23-P1P2.txt")
gamma.r2.parp1.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/PARP1_+_226595744.23-P1P2.txt")
gamma.r2.parp1.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/PARP1_+_226595771.23-P1P2.txt")

gamma.r1.a.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/RNASEH2A_-_12917521.23-P1P2.txt")
gamma.r1.a.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/RNASEH2A_+_12917430.23-P1P2.txt")
gamma.r2.a.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/RNASEH2A_-_12917521.23-P1P2.txt")
gamma.r2.a.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/RNASEH2A_+_12917430.23-P1P2.txt")

gamma.r1.b.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/RNASEH2B_+_51483920.23-P1P2.txt")
gamma.r1.b.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/RNASEH2B_+_51483924.23-P1P2.txt")
gamma.r2.b.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/RNASEH2B_+_51483920.23-P1P2.txt")
gamma.r2.b.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/RNASEH2B_+_51483924.23-P1P2.txt")

gamma.r1.c.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/RNASEH2C_-_65488303.23-P1P2.txt")
gamma.r1.c.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/RNASEH2C_+_65488196.23-P1P2.txt")
gamma.r2.c.1 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/RNASEH2C_-_65488303.23-P1P2.txt")
gamma.r2.c.2 <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/RNASEH2C_+_65488196.23-P1P2.txt")

gamma.r1.all <- read_tsv("../DATA/Interaction_Scores/Compiled/GuideCombination_Scores/guide_combination_interaction_scores_gamma_r1.txt")
gamma.r2.all <- read_tsv("../DATA/Interaction_Scores/Compiled/GuideCombination_Scores/guide_combination_interaction_scores_gamma_r2.txt")
gamma.avg.all <- read_tsv("../DATA/Interaction_Scores/Compiled/GuideCombination_Scores/guide_combination_interaction_scores_gamma_oi_avg.txt")

# Import tau construct scores
tau.avg.all <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_tau_oi_avg.txt")

# Import interaction labels
ints <- read_tsv("../DATA/Interaction_Scores/InteractionCategories/genecombination_interactioncategories_all_interactions.txt")

# Import id map
id.map <- read_tsv("../DATA/guidecombination_id_map.txt")

# Merge in ids
gamma.r1.all <- inner_join(gamma.r1.all, id.map)
gamma.r2.all <- inner_join(gamma.r2.all, id.map)
gamma.avg.all <- inner_join(gamma.avg.all, id.map)

# Combine info from PARP1 models by replicate
gamma.r1.parp1 <- rbind(gamma.r1.parp1.1[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")],
                        gamma.r1.parp1.2[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")])
gamma.r1.parp1$repl <- "R1"

gamma.r2.parp1 <- rbind(gamma.r2.parp1.1[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")],
                        gamma.r2.parp1.2[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")])
gamma.r2.parp1$repl <- "R2"

colnames(gamma.r1.parp1)[4] <- "Gamma"
colnames(gamma.r2.parp1)[4] <- "Gamma"

# Bind replicate info together for plotting
df.parp1 <- rbind(gamma.r2.parp1,
                  gamma.r1.parp1)


# Combine info from RNASEH2A models by replicate
gamma.r1.a <- rbind(gamma.r1.a.1[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")],
                        gamma.r1.a.2[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")])
gamma.r1.a$repl <- "R1"

gamma.r2.a <- rbind(gamma.r2.a.1[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")],
                        gamma.r2.a.2[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")])
gamma.r2.a$repl <- "R2"

colnames(gamma.r1.a)[4] <- "Gamma"
colnames(gamma.r2.a)[4] <- "Gamma"

# Bind replicate info together for plotting
df.a <- rbind(gamma.r2.a,
                  gamma.r1.a)



# Combine info from BRCA2 models by replicate
gamma.r1.b <- rbind(gamma.r1.b.1[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")],
                        gamma.r1.b.2[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")])
gamma.r1.b$repl <- "R1"

gamma.r2.b <- rbind(gamma.r2.b.1[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")],
                        gamma.r2.b.2[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")])
gamma.r2.b$repl <- "R2"

colnames(gamma.r1.b)[4] <- "Gamma"
colnames(gamma.r2.b)[4] <- "Gamma"

# Bind replicate info together for plotting
df.b <- rbind(gamma.r2.b,
                  gamma.r1.b)


# Combine info from BRCA2 models by replicate
gamma.r1.c <- rbind(gamma.r1.c.1[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")],
                        gamma.r1.c.2[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")])
gamma.r1.c$repl <- "R1"

gamma.r2.c <- rbind(gamma.r2.c.1[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")],
                        gamma.r2.c.2[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")])
gamma.r2.c$repl <- "R2"

colnames(gamma.r1.c)[4] <- "Gamma"
colnames(gamma.r2.c)[4] <- "Gamma"

# Bind replicate info together for plotting
df.c <- rbind(gamma.r2.c,
                  gamma.r1.c)

###########################################################

parp1.repl.guide.diffs <- ggplot(df.parp1, aes(x = single, y = Gamma, col = repl)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .6, shape = 16) + 
  facet_wrap(~ query) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  scale_color_manual(values = c("black", "gray70")) +
  xlab("Single guide gamma phenotype") +
  ylab("Combinatorial gamma phenotype") +
  labs(col = "Replicate") +
  geom_text(data = guidelab, aes(label = lab), hjust = 0, size = 3)

c.repl.guide.diffs <- ggplot(df.c, aes(x = single, y = Gamma, col = repl)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .6, shape = 16) + 
  facet_wrap(~ query) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  scale_color_manual(values = c("black", "gray70")) +
  xlab("Single guide gamma phenotype") +
  ylab("Combinatorial gamma phenotype") +
  labs(col = "Replicate") +
  geom_text(data = guidelab, aes(label = lab), hjust = 0, size = 3)