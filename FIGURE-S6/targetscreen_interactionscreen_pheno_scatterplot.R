library(readr)
library(dplyr)
library(ggplot2)

kphenos <- read_tsv("../DATA/AUNIP_targeted/HTS_JL021/aunip_targeted_library_all_phenotypes.txt")
rphenos <- read_tsv("../DATA/AUNIP_targeted/HTS_JL020/aunip_targeted_library_all_phenotypes.txt")

screen.phenos <- read_tsv("../DATA/raw_phenotypes.txt")
id.map <- read_tsv("../DATA/genecombination_id_map.txt")
screen.gis.all.info <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")
screen.gene.phenos <- screen.phenos %>% group_by(GeneCombinationID) %>%
  filter(Identical == FALSE) %>%
  summarise(Gamma.mean = mean(Gamma.Avg),
            Gamma.sd = sd(Gamma.Avg),
            Tau.mean = mean(Tau.Avg),
            Tau.sd = sd(Tau.Avg))
screen.gene.phenos <- inner_join(screen.gene.phenos, id.map)


kphenos.avg <- cbind(rowMeans(kphenos[,c("Gamma.R1", "Gamma.R2", "Gamma.R3")]),
                    rowMeans(kphenos[,c("NIRAP.Tau.R1", "NIRAP.Tau.R2", "NIRAP.Tau.R3")]),
                    rowMeans(kphenos[,c("VELIP1.Tau.R1", "VELIP1.Tau.R2", "VELIP1.Tau.R3")]),
                    rowMeans(kphenos[,c("VELIP2.Tau.R1", "VELIP2.Tau.R2", "VELIP2.Tau.R3")]),
                    rowMeans(kphenos[,c("TALAZOP.Tau.R1", "TALAZOP.Tau.R2", "TALAZOP.Tau.R3")]),
                    rowMeans(kphenos[,c("HU25.Tau.R1", "HU25.Tau.R2", "HU25.Tau.R3")]),
                    rowMeans(kphenos[,c("HU100.Tau.R1", "HU100.Tau.R2", "HU100.Tau.R3")]),
                    rowMeans(kphenos[,c("NIRAP.Rho.R1", "NIRAP.Rho.R2", "NIRAP.Rho.R3")]),
                    rowMeans(kphenos[,c("VELIP1.Rho.R1", "VELIP1.Rho.R2", "VELIP1.Rho.R3")]),
                    rowMeans(kphenos[,c("VELIP2.Rho.R1", "VELIP2.Rho.R2", "VELIP2.Rho.R3")]),
                    rowMeans(kphenos[,c("TALAZOP.Rho.R1", "TALAZOP.Rho.R2", "TALAZOP.Rho.R3")]),
                    rowMeans(kphenos[,c("HU25.Rho.R1", "HU25.Rho.R2", "HU25.Rho.R3")]),
                    rowMeans(kphenos[,c("HU100.Rho.R1", "HU100.Rho.R2", "HU100.Rho.R3")]))
colnames(kphenos.avg) <- c("Gamma", "NIRAP.Tau", "VELIP1.Tau", "VELIP2.Tau", 
                          "TALAZOP.Tau", "HU25.Tau", "HU100.Tau", 
                          "NIRAP.Rho", "VELIP1.Rho", "VELIP2.Rho", 
                          "TALAZOP.Rho", "HU25.Rho", "HU100.Rho")
kphenos.avg <- cbind(kphenos[,1:8], kphenos.avg)
kphenos.avg$Identical <- kphenos.avg$FirstPosition == kphenos.avg$SecondPosition

kgene.phenos <- kphenos.avg %>% group_by(GeneCombinationID) %>%
  filter(Identical == FALSE) %>%
  summarise(Gamma.mean = mean(Gamma),
            Gamma.sd = sd(Gamma),
            NIRAP.Tau.mean = mean(NIRAP.Tau),
            NIRAP.Tau.sd = sd(NIRAP.Tau),
            VELIP1.Tau.mean = mean(VELIP1.Tau),
            VELIP1.Tau.sd = sd(VELIP1.Tau),
            VELIP2.Tau.mean = mean(VELIP2.Tau),
            VELIP2.Tau.sd = sd(VELIP2.Tau),
            TALAZOP.Tau.mean = mean(TALAZOP.Tau),
            TALAZOP.Tau.sd = sd(TALAZOP.Tau),
            HU25.Tau.mean = mean(HU25.Tau),
            HU25.Tau.sd = sd(HU25.Tau),
            HU100.Tau.mean = mean(HU100.Tau),
            HU100.Tau.sd = sd(HU100.Tau),
            NIRAP.Rho.mean = mean(NIRAP.Rho),
            NIRAP.Rho.sd = sd(NIRAP.Rho),
            VELIP1.Rho.mean = mean(VELIP1.Rho),
            VELIP1.Rho.sd = sd(VELIP1.Rho),
            VELIP2.Rho.mean = mean(VELIP2.Rho),
            VELIP2.Rho.sd = sd(VELIP2.Rho),
            TALAZOP.Rho.mean = mean(TALAZOP.Rho),
            TALAZOP.Rho.sd = sd(TALAZOP.Rho),
            HU25.Rho.mean = mean(HU25.Rho),
            HU25.Rho.sd = sd(HU25.Rho),
            HU100.Rho.mean = mean(HU100.Rho),
            HU100.Rho.sd = sd(HU100.Rho),
            N = n())

kgene.phenos <- inner_join(kgene.phenos, unique(kphenos.avg[,c("GeneCombinationID", "FirstGene", "SecondGene")]))
kgene.phenos2 <- kgene.phenos %>% group_by(GeneCombinationID) %>% mutate(sel = min(FirstGene) == FirstGene)
kgene.phenos2 <- kgene.phenos2[kgene.phenos2$sel,]



kcmp.phenos <- kgene.phenos2[,c(1:5,28:30)]
kcmp.phenos$GeneCombinationName <- paste(kcmp.phenos$FirstGene, kcmp.phenos$SecondGene, sep = ":")

kscreen.gene.phenos2 <- screen.gene.phenos[screen.gene.phenos$GeneCombinationName %in% kcmp.phenos$GeneCombinationName,]

kcmp.phenos2 <- inner_join(kscreen.gene.phenos2[,c(6,7,2:5)], kcmp.phenos, 
                          by = c("Gene1" = "FirstGene", "Gene2" = "SecondGene"))

kcmp.phenos2$AUNIP <- grepl("AUNIP", kcmp.phenos2$Gene1)
kcmp.phenos2$NT <- grepl("non-targeting", kcmp.phenos2$Gene1) | grepl("non-targeting", kcmp.phenos2$Gene2)
kcmp.phenos2$OTHER <- kcmp.phenos2$Gene2
kcmp.phenos2$OTHER[kcmp.phenos2$OTHER == "non-targeting"] <- kcmp.phenos2$Gene1[kcmp.phenos2$OTHER == "non-targeting"]

###########

rphenos.avg <- cbind(rowMeans(rphenos[,c("Gamma.R1", "Gamma.R2", "Gamma.R3")]),
                     rowMeans(rphenos[,c("NIRAP.Tau.R1", "NIRAP.Tau.R2", "NIRAP.Tau.R3")]),
                     rowMeans(rphenos[,c("VELIP.Tau.R1", "VELIP.Tau.R2", "VELIP.Tau.R3")]),
                     rowMeans(rphenos[,c("TALAZOP.Tau.R1", "TALAZOP.Tau.R2", "TALAZOP.Tau.R3")]),
                     rowMeans(rphenos[,c("HU100.Tau.R1", "HU100.Tau.R2", "HU100.Tau.R3")]),
                     rowMeans(rphenos[,c("HU400.Tau.R1", "HU400.Tau.R2", "HU400.Tau.R3")]),
                     rowMeans(rphenos[,c("NIRAP.Rho.R1", "NIRAP.Rho.R2", "NIRAP.Rho.R3")]),
                     rowMeans(rphenos[,c("VELIP.Rho.R1", "VELIP.Rho.R2", "VELIP.Rho.R3")]),
                     rowMeans(rphenos[,c("TALAZOP.Rho.R1", "TALAZOP.Rho.R2", "TALAZOP.Rho.R3")]),
                     rowMeans(rphenos[,c("HU100.Rho.R1", "HU100.Rho.R2", "HU100.Rho.R3")]),
                     rowMeans(rphenos[,c("HU400.Rho.R1", "HU400.Rho.R2", "HU400.Rho.R3")]))
colnames(rphenos.avg) <- c("Gamma", "NIRAP.Tau", "VELIP.Tau", 
                           "TALAZOP.Tau", "HU100.Tau", "HU400.Tau", 
                           "NIRAP.Rho", "VELIP.Rho",  
                           "TALAZOP.Rho", "HU100.Rho", "HU400.Rho")
rphenos.avg <- cbind(rphenos[,1:8], rphenos.avg)
rphenos.avg$Identical <- rphenos.avg$FirstPosition == rphenos.avg$SecondPosition

rgene.phenos <- rphenos.avg %>% group_by(GeneCombinationID) %>%
  filter(Identical == FALSE) %>%
  summarise(Gamma.mean = mean(Gamma),
            Gamma.sd = sd(Gamma),
            NIRAP.Tau.mean = mean(NIRAP.Tau),
            NIRAP.Tau.sd = sd(NIRAP.Tau),
            VELIP.Tau.mean = mean(VELIP.Tau),
            VELIP.Tau.sd = sd(VELIP.Tau),
            TALAZOP.Tau.mean = mean(TALAZOP.Tau),
            TALAZOP.Tau.sd = sd(TALAZOP.Tau),
            HU100.Tau.mean = mean(HU100.Tau),
            HU100.Tau.sd = sd(HU100.Tau),
            HU400.Tau.mean = mean(HU400.Tau),
            HU400.Tau.sd = sd(HU400.Tau),
            NIRAP.Rho.mean = mean(NIRAP.Rho),
            NIRAP.Rho.sd = sd(NIRAP.Rho),
            VELIP.Rho.mean = mean(VELIP.Rho),
            VELIP.Rho.sd = sd(VELIP.Rho),
            TALAZOP.Rho.mean = mean(TALAZOP.Rho),
            TALAZOP.Rho.sd = sd(TALAZOP.Rho),
            HU100.Rho.mean = mean(HU100.Rho),
            HU100.Rho.sd = sd(HU100.Rho),
            HU400.Rho.mean = mean(HU400.Rho),
            HU400.Rho.sd = sd(HU400.Rho),
            N = n())

rgene.phenos <- inner_join(rgene.phenos, unique(rphenos.avg[,c("GeneCombinationID", "FirstGene", "SecondGene")]))
rgene.phenos2 <- rgene.phenos %>% group_by(GeneCombinationID) %>% mutate(sel = min(FirstGene) == FirstGene)
rgene.phenos2 <- rgene.phenos2[rgene.phenos2$sel,]



rcmp.phenos <- rgene.phenos2[,c(1:5,24:26)]
rcmp.phenos$GeneCombinationName <- paste(rcmp.phenos$FirstGene, rcmp.phenos$SecondGene, sep = ":")

rscreen.gene.phenos2 <- screen.gene.phenos[screen.gene.phenos$GeneCombinationName %in% rcmp.phenos$GeneCombinationName,]

rcmp.phenos2 <- inner_join(rscreen.gene.phenos2[,c(6,7,2:5)], rcmp.phenos, 
                           by = c("Gene1" = "FirstGene", "Gene2" = "SecondGene"))

rcmp.phenos2$AUNIP <- grepl("AUNIP", rcmp.phenos2$Gene1)
rcmp.phenos2$NT <- grepl("non-targeting", rcmp.phenos2$Gene1) | grepl("non-targeting", rcmp.phenos2$Gene2)
rcmp.phenos2$OTHER <- rcmp.phenos2$Gene2
rcmp.phenos2$OTHER[rcmp.phenos2$OTHER == "non-targeting"] <- rcmp.phenos2$Gene1[rcmp.phenos2$OTHER == "non-targeting"]

#########

plt.k <- ggplot(kcmp.phenos2[kcmp.phenos2$OTHER %in% c("BRCC3", "UIMC1", "FAM175A"),], 
                aes(shape = NT)) +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_abline(alpha = .5) +
  geom_point(aes(Gamma.mean.x, Gamma.mean.y), size = 5, col = "darkred") + 
  geom_point(aes(Tau.mean, NIRAP.Tau.mean), size = 5, col = "darkorchid4") + 
  geom_point(data = kcmp.phenos2[(kcmp.phenos2$OTHER == "AUNIP" & kcmp.phenos2$NT),],
             aes(Gamma.mean.x, Gamma.mean.y), size = 5, stroke = 2, col = "#999999") +
  geom_point(data = kcmp.phenos2[(kcmp.phenos2$OTHER == "AUNIP" & kcmp.phenos2$NT),],
             aes(Tau.mean, NIRAP.Tau.mean), size = 5, stroke = 2, col = "black") +             
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        aspect.ratio = 2/1) +
  scale_shape_manual(name = "Other Gene",
                     labels = c("AUNIP", "NT"),
                     values = c(16, 3)) +
  # scale_color_manual(name = "Gene",
  #                    values = c("black", "#004949", "#ff6db6", "#490092", "#b66dff", "#b6dbff", "#920000", "#db6d00")) +
  # guides(color = guide_legend(order = 1),
  #        shape = guide_legend(order = 2)) +
  coord_fixed() +
  xlab("Interaction Screen Phenotype") + 
  ylab("Target Screen K562 Phenotype") +
  scale_y_continuous(limits = c(-.46, .16),
                     breaks = seq(-.4, .1, .1),
                     labels = c(-.4, -.3, -.2, -.1, 0, .1)) +
  scale_x_continuous(limits = c(-.22, .01),
                     breaks = seq(-.2, 0, .1),
                     labels = c(-.2, -.1, 0))

plt.r <- ggplot(rcmp.phenos2[rcmp.phenos2$OTHER %in% c("BRCC3", "UIMC1", "FAM175A"),], 
                aes(shape = NT)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_abline() +
  geom_point(aes(Gamma.mean.x, Gamma.mean.y), size = 5, col = "darkred") + 
  geom_point(aes(Tau.mean, NIRAP.Tau.mean), size = 5, col = "darkorchid4") + 
  geom_point(data = rcmp.phenos2[(rcmp.phenos2$OTHER == "AUNIP" & rcmp.phenos2$NT),],
             aes(Gamma.mean.x, Gamma.mean.y), size = 5, stroke = 2, col = "#999999") +
  geom_point(data = rcmp.phenos2[(rcmp.phenos2$OTHER == "AUNIP" & rcmp.phenos2$NT),],
             aes(Tau.mean, NIRAP.Tau.mean), size = 5, stroke = 2, col = "black") +             
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none",
        aspect.ratio = 2/1) +
  scale_shape_manual(name = "Other Gene",
                     labels = c("AUNIP", "NT"),
                     values = c(16, 3)) +
  # scale_color_manual(name = "Gene",
  #                    values = c("black", "#004949", "#ff6db6", "#490092", "#b66dff", "#b6dbff", "#920000", "#db6d00")) +
  # guides(color = guide_legend(order = 1),
  #        shape = guide_legend(order = 2)) +
  xlab("Interaction Screen Phenotype") + 
  ylab("Target Screen RPE1 Phenotype") +
  scale_y_continuous(limits = c(-.46, .16),
                     breaks = seq(-.4, .1, .1),
                     labels = c(-.4, -.3, -.2, -.1, 0, .1)) +
  scale_x_continuous(limits = c(-.22, .01),
                     breaks = seq(-.2, 0, .1),
                     labels = c(-.2, -.1, 0))

plt <- plot_grid(plt.k, plt.r)

ggsave("../FIGURES/FIGURE-6/target_interaction_screen_phenotype_comparison.svg", plt,
       height = 6, width = 6, device = "svg")
