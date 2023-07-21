library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)

# THESE ARE THE RPE1 DATA BUT I'M LEAVING EVERYTHING K CAUSE LAZY
kphenos <- read_tsv("../DATA/AUNIP_targeted/HTS_JL020/aunip_targeted_library_all_phenotypes.txt")

kphenos2 <- kphenos %>% 
  filter(GeneCombinationID %in% c("tlgc2", "tlgc4", "tlgc7")) %>% 
  select(GuideCombinationID, 
         Gamma.R1:Gamma.R3, 
         NIRAP.Rho.R1:HU100.Rho.R3) %>%
  group_by(GuideCombinationID) %>% 
  summarise(across(everything(), list(mean = mean)))

kphenos2 <- inner_join(kphenos2, unique(kphenos[,c("GuideCombinationID", "GeneCombinationID")]))
kphenos2 <- kphenos2[,c(1,17,2:16)]
kphenos3 <- kphenos2[,2:17] %>% 
  group_by(GeneCombinationID) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))


k.melt <- melt(kphenos3)
k.melt$Phenotype <- gsub("\\..*", "", k.melt$variable)
k.melt$Phenotype[k.melt$Phenotype == "Gamma"] <- "GAMMA"
k.melt$Phenotype[k.melt$Phenotype == "HU100"] <- "HU"

k.melt$Gene <- "BRCC3"
k.melt$Gene[k.melt$GeneCombinationID == "tlgc4"] <- "FAM175A"
k.melt$Gene[k.melt$GeneCombinationID == "tlgc7"] <- "UIMC1"

k.melt$measure <- "mean"
k.melt$measure[grepl("sd", k.melt$variable)] <- "sd"

k.melt$rep <- gsub(".*\\.", "", k.melt$variable)
k.melt$rep <- gsub("_.*", "", k.melt$rep)

k.wide <- reshape(k.melt[,c(1,3:7)],
                  idvar = c("GeneCombinationID", "Phenotype", "Gene", "rep"),
                  timevar = "measure",
                  direction = "wide")

k.avg <- k.wide %>% 
  group_by(Gene, Phenotype) %>%
  summarise(avg = mean(value.mean))
k.avg$variable <- "OBSERVED"
colnames(k.avg)[3] <- "value"
k.avg <- k.avg[,c("Gene", "Phenotype", "variable", "value")]


ks.phenos <- kphenos %>% 
  filter(GeneCombinationID %in% c("tlgc8", "tlgc9", "tlgc11", "tlgc14")) %>% 
  select(GuideCombinationID, 
         Gamma.R1:Gamma.R3, 
         NIRAP.Rho.R1:HU100.Rho.R3) %>%
  group_by(GuideCombinationID) %>% 
  summarise(across(everything(), list(mean = mean)))

ks.phenos <- inner_join(ks.phenos, unique(kphenos[,c("GuideCombinationID", "GeneCombinationID")]))
ks.phenos <- ks.phenos[,c(1,17,2:16)]
ks.phenos2 <- ks.phenos[,2:17] %>% 
  group_by(GeneCombinationID) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))


ks.melt <- melt(ks.phenos2)
ks.melt$Phenotype <- gsub("\\..*", "", ks.melt$variable)
ks.melt$Phenotype[ks.melt$Phenotype == "Gamma"] <- "GAMMA"
ks.melt$Phenotype[ks.melt$Phenotype == "HU100"] <- "HU"

ks.melt$Gene <- "BRCC3"
ks.melt$Gene[ks.melt$GeneCombinationID == "tlgc8"] <- "AUNIP"
ks.melt$Gene[ks.melt$GeneCombinationID == "tlgc11"] <- "FAM175A"
ks.melt$Gene[ks.melt$GeneCombinationID == "tlgc14"] <- "UIMC1"

ks.melt$measure <- "mean"
ks.melt$measure[grepl("sd", ks.melt$variable)] <- "sd"

ks.melt$rep <- gsub(".*\\.", "", ks.melt$variable)
ks.melt$rep <- gsub("_.*", "", ks.melt$rep)

ks.wide <- reshape(ks.melt[,c(1,3:7)],
                   idvar = c("GeneCombinationID", "Phenotype", "Gene", "rep"),
                   timevar = "measure",
                   direction = "wide")

ks.avg <- ks.wide %>% 
  group_by(Gene, Phenotype) %>%
  summarise(avg = mean(value.mean))

ks.avg2 <- ks.avg[ks.avg$Gene != "AUNIP",]
ks.avg2$AUNIP <- rep(ks.avg$avg[ks.avg$Gene == "AUNIP"], 3)
ks.avg2$ADDITIVE <- ks.avg2$avg + ks.avg2$AUNIP
colnames(ks.avg2)[3] <- "SINGLE"

ks.avg.melt <- melt(data = ks.avg2)

#### BARS FOR PLOT
k.bars <- rbind(k.avg, ks.avg.melt)
k.bars$Gene <- factor(k.bars$Gene, levels = c("BRCC3", "FAM175A", "UIMC1"))
k.bars$Phenotype <- factor(k.bars$Phenotype, levels = c("GAMMA", "VELIP", "NIRAP", "TALAZOP", "HU"))
k.bars$variable <- factor(k.bars$variable, levels = c("SINGLE", "AUNIP", "ADDITIVE", "OBSERVED"))
k.bars$Int <- interaction(k.bars$Gene, k.bars$variable)
k.bars$Int <- as.character(k.bars$Int)
k.bars$Int[k.bars$variable == "AUNIP"] <- "AUNIP"
k.bars$Int[k.bars$variable == "ADDITIVE"] <- "ADDITIVE"
k.bars$Int[k.bars$variable == "SINGLE"] <- gsub(".SINGLE", "", k.bars$Int[k.bars$variable == "SINGLE"])
k.bars$Int[k.bars$variable == "OBSERVED"] <- gsub(".OBSERVED", "", k.bars$Int[k.bars$variable == "OBSERVED"])
k.bars$Int[k.bars$variable == "OBSERVED"] <- paste0("AUNIP:", k.bars$Int[k.bars$variable == "OBSERVED"])
k.bars$Int <- factor(k.bars$Int, levels = c("BRCC3", "FAM175A", "UIMC1",
                                            "AUNIP", "ADDITIVE",
                                            "AUNIP:BRCC3", "AUNIP:FAM175A", "AUNIP:UIMC1"))

#### PTS WITH ERRORS FOR PLOT
k.wide2 <- k.wide[,c("Gene", "Phenotype", "rep", "value.mean", "value.sd")]
k.wide2$variable <- "OBSERVED"
ks.wide2 <- ks.wide[,c("Gene", "Phenotype", "rep", "value.mean", "value.sd")]

ks.wide3 <- ks.wide2[ks.wide2$Gene != "AUNIP",]
aunip.wide <- ks.wide2[ks.wide2$Gene == "AUNIP",]

ks.wide4 <- left_join(ks.wide3, aunip.wide[,2:5],
                      by = c("Phenotype", "rep"))
colnames(ks.wide4)[4:7] <- c("SINGLE.mean", "SINGLE.sd", "AUNIP.mean", "AUNIP.sd")

ks.wide.melt <- melt(ks.wide4)
ks.wide.melt$measure <- "mean"
ks.wide.melt$measure[grepl("sd",ks.wide.melt$variable)] <- "sd"
ks.wide.melt$variable <- gsub("\\..*", "", ks.wide.melt$variable)
ks.wide.melt.wide <- reshape(ks.wide.melt, 
                             idvar= c("Gene", "Phenotype", "rep", "variable"),
                             timevar= "measure",
                             direction = "wide")
ks.add <- data.frame(Gene = rep(c("BRCC3", "FAM175A", "UIMC1"), each = 5),
                     Phenotype = rep(c("GAMMA", "NIRAP", "VELIP", "TALAZOP", "HU"), 3),
                     variable = "ADDITIVE",
                     rep = "R1",
                     value.mean = NA,
                     value.sd = NA)

k.pts <- rbind(k.wide2[,c("Gene", "Phenotype", "variable", "rep", "value.mean", "value.sd")],
               ks.wide.melt.wide[,c("Gene", "Phenotype", "variable", "rep", "value.mean", "value.sd")],
               ks.add)
k.pts$Int <- interaction(k.pts$Gene, k.pts$variable)
k.pts$Int <- as.character(k.pts$Int)
k.pts$Int[k.pts$variable == "AUNIP"] <- "AUNIP2"
k.pts$Int[k.pts$variable == "ADDITIVE"] <- "ADDITIVE"
k.pts$Int[k.pts$variable == "SINGLE"] <- paste0(gsub(".SINGLE", "", k.pts$Int[k.pts$variable == "SINGLE"]),2)
k.pts$Int[k.pts$variable == "OBSERVED"] <- gsub(".OBSERVED", "", k.pts$Int[k.pts$variable == "OBSERVED"])
k.pts$Int[k.pts$variable == "OBSERVED"] <- paste0("AUNIP:", k.pts$Int[k.pts$variable == "OBSERVED"], "2")
k.pts$Int <- factor(k.pts$Int, levels = c("BRCC32", "FAM175A2", "UIMC12",
                                          "AUNIP2", "ADDITIVE",
                                          "AUNIP:BRCC32", "AUNIP:FAM175A2", "AUNIP:UIMC12"))

k.range <- k.pts %>%
  group_by(Gene, Phenotype, variable) %>%
  summarise(min = min(value.mean),
            max = max(value.mean))
k.range <- inner_join(k.range, k.bars[,c("Gene","Phenotype", "variable", "Int")],
                      by = c("Gene", "Phenotype", "variable"))

k.bars$Phenotype <- factor(k.bars$Phenotype, levels = c("GAMMA", "VELIP", "NIRAP", "TALAZOP", "HU"))
k.pts$Phenotype <- factor(k.pts$Phenotype, levels = c("GAMMA", "VELIP", "NIRAP", "TALAZOP", "HU"))
k.range$Phenotype <- factor(k.range$Phenotype, levels = c("GAMMA", "VELIP", "NIRAP", "TALAZOP", "HU"))


p1 <- ggplot(k.bars[k.bars$Phenotype == "GAMMA",], aes(Gene)) + 
  geom_hline(yintercept = .25, linetype = "dotted") +
  geom_hline(yintercept = -.25, linetype = "dotted") +
  geom_hline(yintercept = -.5, linetype = "dotted") +
  geom_hline(yintercept = -.75, linetype = "dotted") +
  geom_hline(yintercept = -1, linetype = "dotted") +
  geom_hline(yintercept = -1.25, linetype = "dotted") +
  geom_hline(yintercept = -1.5, linetype = "dotted") +
  geom_hline(yintercept = -1.75, linetype = "dotted") +
  facet_wrap(~ Phenotype, 
             nrow = 1,
             strip.position = "bottom") + 
  geom_col(aes(y = value, fill = Int),
           position = position_dodge()) +
  geom_linerange(data = k.range[k.range$Phenotype == "GAMMA",],
                 aes(ymin = min, ymax = max, group = Int),
                 position = position_dodge(width = .9)) +
  geom_point(data = k.pts[k.pts$Phenotype == "GAMMA",], 
             aes(y = value.mean, fill = Int),
             shape = 21, color = "black", size = 3, 
             position = position_dodge(width = .9)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none") +
  xlab("") + ylab("Gamma Phenotype") +
  scale_fill_manual(values = c("#333333", "#88ccee", 
                               "#997700","#b59322",  "#225522",    
                               "#3a763a", "#994455",  "#b56071",  "#7AB7D6",
                               "#eecc66", "#d2b044" ,  "#69b869",  
                               "#519751",  "#ee99aa",  "#D27D8E")) 

p2 <- ggplot(k.bars[k.bars$Phenotype != "GAMMA",], aes(Gene)) + 
  geom_hline(yintercept = .25, linetype = "dotted") +
  geom_hline(yintercept = -.25, linetype = "dotted") +
  geom_hline(yintercept = -.5, linetype = "dotted") +
  geom_hline(yintercept = -.75, linetype = "dotted") +
  geom_hline(yintercept = -1, linetype = "dotted") +
  geom_hline(yintercept = -1.25, linetype = "dotted") +
  geom_hline(yintercept = -1.5, linetype = "dotted") +
  geom_hline(yintercept = -1.75, linetype = "dotted") +
  facet_wrap(~ Phenotype, 
             nrow = 1,
             strip.position = "bottom") + 
  geom_col(aes(y = value, fill = Int),
           position = position_dodge()) +
  geom_linerange(data = k.range[k.range$Phenotype != "GAMMA",],
                 aes(ymin = min, ymax = max, group = Int),
                 position = position_dodge(width = .9)) +
  geom_point(data = k.pts[k.pts$Phenotype != "GAMMA",], 
             aes(y = value.mean, fill = Int),
             shape = 21, color = "black", size = 3,
             position = position_dodge(width = .9)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none") +
  xlab("") + ylab("Rho Phenotype") +
  scale_fill_manual(values = c("#333333", "#88ccee", 
                               "#997700","#b59322",  "#225522",    
                               "#3a763a", "#994455",  "#b56071",  "#7AB7D6",
                               "#eecc66", "#d2b044" ,  "#69b869",  
                               "#519751",  "#ee99aa",  "#D27D8E"))   

p3 <- plot_grid(p1, p2, rel_widths = c(1,3))

ggsave("../FIGURES/FIGURE-6/aunip_screen_rpe1_rho_sidebyside.svg", p3,
       height = 5, width = 7)
