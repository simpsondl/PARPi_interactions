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


nup <- "NUP62_+_50432687.23-P1P2"
gamma.r1.nup <- read_tsv("../DATA/Interaction_Scores/Gamma.R1/NUP62_+_50432687.23-P1P2.txt")
gamma.r2.nup <- read_tsv("../DATA/Interaction_Scores/Gamma.R2/NUP62_+_50432687.23-P1P2.txt")
gamma.avg.nup <- read_tsv("../DATA/Interaction_Scores/Gamma.Avg/NUP62_+_50432687.23-P1P2.txt")
gamma.oiavg.nup <- read_tsv("../DATA/Interaction_Scores/Gamma.OI.Avg/NUP62_+_50432687.23-P1P2.txt")

model_coeffs_r1 <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_gamma_r1.txt")
model_coeffs_r1 <- model_coeffs_r1[model_coeffs_r1$sgRNA.ID == nup,]
model_coeffs_r2 <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_gamma_r2.txt")
model_coeffs_r2 <- model_coeffs_r2[model_coeffs_r2$sgRNA.ID == nup,]
model_coeffs_avg <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_gamma_avg.txt")
model_coeffs_avg <- model_coeffs_avg[model_coeffs_avg$sgRNA.ID == nup,]
model_coeffs_oiavg <- read_tsv("../DATA/Interaction_Scores/Compiled/Model_Results/model_estimates_gamma_oi_avg.txt")
model_coeffs_oiavg <- model_coeffs_oiavg[model_coeffs_oiavg$sgRNA.ID == nup,]


gamma.r1.nup <- gamma.r1.nup[,c("sgRNA.id", "query", "single", "Gamma.R1", "GI.z")]
gamma.r1.nup$repl <- "R1"

gamma.r2.nup <- gamma.r2.nup[,c("sgRNA.id", "query", "single", "Gamma.R2", "GI.z")]
gamma.r2.nup$repl <- "R2"

gamma.avg.nup <- gamma.avg.nup[,c("sgRNA.id", "query", "single", "Gamma.Avg", "GI.z")]
gamma.avg.nup$repl <- "AVG"

gamma.oiavg.nup <- gamma.oiavg.nup[,c("sgRNA.id", "query", "single", "Gamma.OI.Avg", "GI.z")]
gamma.oiavg.nup$repl <- "OI.AVG"

colnames(gamma.r1.nup)[4] <- "Gamma"
colnames(gamma.r2.nup)[4] <- "Gamma"
colnames(gamma.avg.nup)[4] <- "Gamma"
colnames(gamma.oiavg.nup)[4] <- "Gamma"

# Bind replicate info together for plotting
df.nup <- rbind(gamma.r2.nup,
                gamma.r1.nup,
                gamma.avg.nup,
                gamma.oiavg.nup)

model_line <- data.frame(single = seq(-.37, .07, .01))
model_line$pred.r1 <- model_coeffs_r1$estimate[1] + 
  model_coeffs_r1$estimate[2] * model_line$single + 
  model_coeffs_r1$estimate[3] * model_line$single^2
model_line$pred.r2 <- model_coeffs_r2$estimate[1] + 
  model_coeffs_r2$estimate[2] * model_line$single + 
  model_coeffs_r2$estimate[3] * model_line$single^2
model_line$pred.avg <- model_coeffs_oiavg$estimate[1] + 
  model_coeffs_oiavg$estimate[2] * model_line$single + 
  model_coeffs_oiavg$estimate[3] * model_line$single^2

df.nup$control <- grepl("non-targeting", df.nup$query) | grepl("non-targeting", df.nup$sgRNA.id)

df.nup$gene <- gsub("_.*", "", df.nup$sgRNA.id)
df.nup.sum <- df.nup %>% 
  group_by(gene, repl) %>% 
  summarise(GI.z = mean(GI.z), control = unique(control)) %>% 
  filter(repl == "OI.AVG")

df.nup$repl <- factor(df.nup$repl, levels = c("R1", "R2", "AVG", "OI.AVG"))

nup.model.sig <- df.nup.sum$gene[abs(df.nup.sum$GI.z) >= 1.524]
df.nup$sig <- df.nup$gene %in% nup.model.sig

# nup.segs <- data.frame(guide = unique(df.nup$sgRNA.id[df.nup$sig]),
#                        xmin = NA,
#                        xmax = NA,
#                        ymin = NA,
#                        ymax = NA)
# for(i in 1:24){
#   nup.segs[i,2:5] <- c(df.nup$single[df.nup$sgRNA.id == nup.segs$guide[i] & df.nup$repl == "R1"],
#                        df.nup$single[df.nup$sgRNA.id == nup.segs$guide[i] & df.nup$repl == "R2"],
#                        df.nup$Gamma[df.nup$sgRNA.id == nup.segs$guide[i] & df.nup$repl == "R1"],
#                        df.nup$Gamma[df.nup$sgRNA.id == nup.segs$guide[i] & df.nup$repl == "R2"])
# }

nup.models <- ggplot(df.nup[df.nup$repl != "OI.AVG",]) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_line(data = model_line, aes(single, pred.avg), size = 1.25) +
  geom_point(data = df.nup[df.nup$repl == "AVG",],
             aes(x = single, y = Gamma),
             alpha = .7, shape = 16, size = 1.5, col = "black") +
  # geom_point(data = df.nup[df.nup$sig & df.nup$repl == "AVG",], 
  #            aes(x = single, y = Gamma),
  #            shape = 16, size = 3, col = "red") + 
  geom_boxplot(data = df.nup[df.nup$control & df.nup$repl == "OI.AVG",], aes(x = 0, y = Gamma),
               width = .017, alpha = .5, size = 1.5, fill = "#33bbee", color = "#33bbee") +
  #geom_segment(data = nup.segs, aes(x = xmin, xend = xmax, y = ymin, yend = ymax), alpha = .4) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.background = element_blank()) +
  xlab("Single guide gamma phenotype") +
  ylab("Combinatorial gamma phenotype with NUP62") +
  labs(col = "Replicate") 

ggsave("../FIGURES/TechnicalNote/nup62_model.pdf", nup.models,
       device = "pdf", height = 6, width = 6)
ggsave("../FIGURES/TechnicalNote/nup62_model.svg", nup.models,
       device = "svg", height = 6, width = 6)

nup.cons <- all.gamma[all.gamma$sgRNA.id == nup | all.gamma$query == nup,]
nup.cons$nonnup <- nup.cons$query
nup.cons$nonnup[nup.cons$nonnup == nup] <- nup.cons$sgRNA.id[nup.cons$nonnup == nup]
nup.cons$isquery <- nup.cons$query == nup
nup.cons$Avg <- rowMeans(nup.cons[,c("Gamma.R1", "Gamma.R2")])


nup.orient <- dcast(nup.cons[,c("isquery", "nonnup", "GI.z")], nonnup ~ isquery)
nup.orient <- nup.orient[!is.na(nup.orient$`TRUE`) & !is.na(nup.orient$`FALSE`),]
nup.orient2 <- inner_join(nup.orient, phenos[,c("sgRNA.ID", "Gamma.Avg")], by = c("nonnup" = "sgRNA.ID"))
nup.orient2$Diff <- nup.orient2$Gamma.Avg - phenos$Gamma.Avg[phenos$sgRNA.ID == "NUP62_+_50432687.23-P1P2"]

nup.fit <- lm(`FALSE` ~ `TRUE`, data = nup.orient2)
nup.fit2 <- lm(`FALSE` ~ `TRUE`, data = nup.orient2[abs(nup.orient2$Gamma.Avg - .1) <= .2,])

nup.orient.plt <- ggplot(nup.orient2, aes(`TRUE`, `FALSE`)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  #geom_point(aes(col = Diff)) +
  geom_point(alpha = .6) +
  theme_bw() + 
  # scale_color_gradientn(colors = c("#543005", "#8C510A", "#BF812D", "#DFC27D",
  #                                  "#bbbbbb", 
  #                                  "#80CDC1", "#35978F", "#01665E", "#003C30"),
  #                       values = c(0, 1/9*1.4, 2/9+ .4*1/9, 3/9 + .4*1/9, 4/9 + .4*1/9,
  #                                  5/9+ .2*1/9,
  #                                  6/9 + .2*1/9, 7/9, 8/9, 1)) +
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

ggsave("../FIGURES/TechnicalNote/nup62_orientation.pdf", nup.orient.plt,
       device = "pdf", height = 6, width = 3)
ggsave("../FIGURES/TechnicalNote/nup62_orientation.svg", nup.orient.plt,
       device = "svg", height = 6, width = 3)

######### Not sure where this is from right now, it's 1:30am and Julia has been crying for 2 hours straight
all.gamma$Avg <- rowMeans(all.gamma[,c("Gamma.R1", "Gamma.R2")])

orient.cor <- data.frame(guide = unique(all.gamma$query), cor = NA )
for(i in orient.cor$guide){
  tmp <- all.gamma[all.gamma$query == i | all.gamma$sgRNA.id == i,]
  tmp$non <- tmp$sgRNA.id
  tmp$non[tmp$non == i] <- tmp$query[tmp$non == i]
  tmp$isquery <- tmp$query == i
  tmp2 <- dcast(tmp[,c("non", "isquery", "Avg")], non ~ isquery)
  orient.cor$cor[orient.cor$guide == i] <- cor(tmp2$`FALSE`, tmp2$`TRUE`, use = "complete.obs")
}

cor.plt <- ggplot(orient.cor, aes(cor)) +
  geom_density(fill = "#bbbbbb", alpha = .5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Correlation between Oriented Phenotypes") + ylab("") +
  annotate("segment", x = 0.18351505, xend = 0.18351505, yend = .25, y = .75,
           arrow = arrow(), size = 1.5) +
  annotate("text", x = 0.18351505, y = .85, label = "NUP62")

ggsave("../FIGURES/TechnicalNote/phenotype_correlations.pdf", cor.plt,
       device = "pdf", height = 6, width = 6, dpi = 300)

nup.gcs <- gamma.gcs[grepl("NUP62", gamma.gcs$GeneCombinationName),]
nup.gcs <- nup.gcs[nup.gcs$N != 80,]

gc.plt <- ggplot(nup.gcs[], aes(as.factor(N), GI.Avg)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_hline(yintercept = nup.fit$coefficients[[1]], alpha = .5, linetype = "dashed") +
  geom_violin() +
  geom_quasirandom(alpha = .5, size = .9) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Number of Supporting Interactions") +
  ylab("Interaction Score")

ggsave("../FIGURES/TechnicalNote/nup_supporting.pdf", gc.plt,
       device = "pdf", height = 6, width = 6, dpi = 300)

gc.var <- gamma.gcs %>% group_by(N) %>% filter(N <= 8) %>% summarise(mean = mean(GI.Avg), var = var(GI.Avg))
gamma.gcs2 <- inner_join(gamma.gcs, gc.var)
gc.counts <- as.data.frame(table(gamma.gcs2$N))

supp.plt <- ggplot(gamma.gcs2, aes(as.factor(N), GI.Avg)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_hline(yintercept = -1.524, alpha = .5, linetype = "dashed") +
  geom_hline(yintercept = 1.524, alpha = .5, linetype = "dashed") +
  geom_violin(aes(fill = var)) +
  geom_text(data = gc.counts, aes(x = as.factor(Var1), y = -16, label = Freq), size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Number of Supporting Interactions") +
  ylab("Interaction Score") + 
  scale_fill_viridis_c(name = "Variance")

ggsave("../FIGURES/TechnicalNote/all_supporting.pdf", supp.plt, device = "pdf",
       height = 6, width = 6, dpi = 300)
