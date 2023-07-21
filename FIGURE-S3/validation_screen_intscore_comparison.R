library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(outliers)
library(effsize)

processflow <- function(fp, col1 = "dodgerblue", col2 = "darkorange1", colexp = "#bbbbbb", colobs = "black") {
  flow <- read_tsv(fp)
  colnames(flow)[1] <- "population"
  flow$replicate <- rep(paste0("R",1:3), each = 4)
  flow <- flow[,c(1,ncol(flow), 2:(ncol(flow)-1))]
  # normalize to day 6 or first timepoint
  flow[,3:(ncol(flow))] <- flow[,3:(ncol(flow))]/as.vector(flow[,3])
  # normalize to -/- population
  flow[1:4,3:ncol(flow)] <- sweep(as.data.frame(flow[1:4,3:ncol(flow)]),
                                  2,
                                  unlist(flow[flow$population == "-/-" & flow$replicate == "R1", 3:ncol(flow)]),
                                  "/")
  flow[5:8,3:ncol(flow)] <- sweep(as.data.frame(flow[5:8,3:ncol(flow)]),
                                  2,
                                  unlist(flow[flow$population == "-/-" & flow$replicate == "R2", 3:ncol(flow)]),
                                  "/")
  flow[9:12,3:ncol(flow)] <- sweep(as.data.frame(flow[9:12,3:ncol(flow)]),
                                   2,
                                   unlist(flow[flow$population == "-/-" & flow$replicate == "R3", 3:ncol(flow)]),
                                   "/")
  # convert to fold-change
  flow[,3:ncol(flow)] <- log2(flow[,3:ncol(flow)])
  # add in additive expected scores
  flow[13,] <- c("Exp", "R1", flow[flow$population == "BFP+" & flow$replicate == "R1", 3:ncol(flow)] +
                   flow[flow$population == "GFP+" & flow$replicate == "R1", 3:ncol(flow)])
  flow[14,] <- c("Exp", "R2", flow[flow$population == "BFP+" & flow$replicate == "R2", 3:ncol(flow)] +
                   flow[flow$population == "GFP+" & flow$replicate == "R2", 3:ncol(flow)])
  flow[15,] <- c("Exp", "R3", flow[flow$population == "BFP+" & flow$replicate == "R3", 3:ncol(flow)] +
                   flow[flow$population == "GFP+" & flow$replicate == "R3", 3:ncol(flow)])
  # set levels for legend
  flow$population <- factor(flow$population, levels = c("BFP+", "GFP+", "Exp", "BFP+GFP+", "-/-"))
  # rearrange for plotting
  flow.m <- melt(flow)
  flow.m$day <- gsub("day", "", flow.m$variable)
  flow.m$day <- as.numeric(flow.m$day)
  # get min/max per group for bars
  flow.error <- flow.m %>% group_by(population, day) %>% summarise(min = min(value), max = max(value))
  flow.error.m <- melt(flow.error, id.vars = c("population", "day"))
  
  # make plot
  plt <- ggplot(flow.m[flow.m$population != "-/-",]) +
    geom_point(aes(day, value, col = population)) +
    geom_line(aes(day, value, 
                  col = population, 
                  group = interaction(population, replicate))) +
    geom_ribbon(data = flow.error[flow.error$population != "-/-",], 
                aes(x = day, ymin = min, ymax = max, fill = population),
                alpha = .5) +
    stat_summary(aes(day, value, col = population), geom = "line", fun = mean, size = 2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.length = unit(.2, "cm")) +
    scale_x_continuous(breaks = unique(flow.m$day)) +
    scale_color_manual(values = c(col1, col2, colexp, colobs)) +
    scale_fill_manual(values = c(col1, col2, colexp, colobs)) +
    ylab("log2 FC")
  
  
  # Test t-test assumptions
  warnings <- NA
  exp <- unlist(flow[flow$population == "Exp",ncol(flow)])
  obs <- unlist(flow[flow$population == "BFP+GFP+",ncol(flow)])
  #check for outliers
  grubbs.exp.max <- grubbs.test(exp)$p.value
  grubbs.exp.min <- grubbs.test(exp, opposite = TRUE)$p.value
  grubbs.obs.max <- grubbs.test(obs)$p.value
  grubbs.obs.min <- grubbs.test(obs, opposite = TRUE)$p.value
  if(grubbs.exp.max <= .05 | grubbs.exp.min <= .05 | grubbs.obs.max <= .05 | grubbs.obs.min <= .05){
    warnings <- "Outlier in dataset"
  }
  # check for normality
  exp.norm <- shapiro.test(exp)$p.value
  obs.norm <- shapiro.test(obs)$p.value
  datnorm <- TRUE
  if(exp.norm <= .05 | obs.norm <= .05){
    warnings <- c(warnings, "Data not normally distributed")
    datnorm <- FALSE
  }
  # Check for equality of variances
  exp.var <- var(exp)
  obs.var <- var(obs)
  vareq <- TRUE
  if(var.test(exp, obs)$p.value <= .05){
    vareq <- FALSE
    warnings <- c(warnings, "Unequal variances")
  }
  ttest <- t.test(exp,
                  obs,
                  paired = TRUE, 
                  var.equal = vareq)$p.value
  eff <- cohen.d(obs, 
                 exp,
                 paired = TRUE,
                 pooled = vareq)$estimate
  
  return(list(warnings, ttest, eff, plt))
}

all.val <- read_tsv("../DATA/screendata_for_validationexperiments.txt")
single.pheno <- read_tsv("../DATA/single_sgRNA_phenotypes.txt")

ints <- read_tsv("../DATA/Interaction_Scores/Compiled/Construct_Scores/all_interaction_scores_gamma_oi_avg.txt")
ints2 <- ints[ints$GuideCombinationID %in% all.val$GuideCombinationID,] %>%
  select(GuideCombinationID, sgRNA.id:single, Gamma.OI.Avg:GI.z)

ints2 <- inner_join(ints2, single.pheno[,c("sgRNA.ID", "Gamma.Avg")],
                    by = c("query" = "sgRNA.ID"))
ints2$factor <- ints2$GI.z/ints2$GI
ints2$Additive <- ints2$single + ints2$Gamma.Avg
ints2$GI.add <- ints2$Gamma.OI.Avg - ints2$Additive
ints2$GI.add.z <- ints2$GI.add*ints2$factor
ints2$Gene1 <- gsub("_.*", "", ints2$sgRNA.id)
ints2$Gene2 <- gsub("_.*", "", ints2$query)
ints2$GeneA <- apply(ints2[,c("Gene1", "Gene2")], 1, min)
ints2$GeneB <- apply(ints2[,c("Gene1", "Gene2")], 1, max)
ints2$int <- paste(ints2$GeneA, ints2$GeneB, sep = ":")

exp.parp1_1.c_1.k562 <- file.path("../DATA/ValidationExperiments/parp1-1-bfp_rnaseh2c-1-gfp_k562.txt")
exp.parp1_2.c_2.k562 <- file.path("../DATA/ValidationExperiments/parp1-2-gfp_rnaseh2c-2-bfp_k562.txt")
exp.parp2_1.c_1.k562 <- file.path("../DATA/ValidationExperiments/parp2-1-bfp_rnaseh2c-1-gfp_k562.txt")
exp.parp2_2.c_2.k562 <- file.path("../DATA/ValidationExperiments/parp2-2-gfp_rnaseh2c-2-bfp_k562.txt")

exp.blm_1.sp_2.k562 <- file.path("../DATA/ValidationExperiments/blm-1-bfp_spidr-2-gfp_k562.txt")
exp.blm_2.sp_1.k562 <- file.path("../DATA/ValidationExperiments/blm-2-gfp_spidr-1-bfp_k562.txt")

exp.mnd_2.hop2_2.k562 <- file.path("../DATA/ValidationExperiments/psmc3ip-2-gfp_mnd1-2-bfp_k562.txt")
exp.mnd_2.rad_2.k562 <- file.path("../DATA/ValidationExperiments/rad54l-2-gfp_mnd1-2-bfp_k562.txt")
exp.hop2_2.rad_2.k562 <- file.path("../DATA/ValidationExperiments/rad54l-2-gfp_psmc3ip-2-bfp_k562.txt")

flow.mnd_2.rad_2.k562 <- processflow(exp.mnd_2.rad_2.k562)
flow.hop2_2.rad_2.k562 <- processflow(exp.hop2_2.rad_2.k562)
flow.mnd_2.hop2_2.k562 <- processflow(exp.mnd_2.hop2_2.k562)

flow.blm_2.sp_1.k562 <- processflow(exp.blm_2.sp_1.k562)
flow.blm_1.sp_2.k562 <- processflow(exp.blm_1.sp_2.k562)

flow.parp1_1.c_1.k562 <- processflow(exp.parp1_1.c_1.k562)
flow.parp1_2.c_2.k562 <- processflow(exp.parp1_2.c_2.k562)
flow.parp2_1.c_1.k562 <- processflow(exp.parp2_1.c_1.k562)
flow.parp2_2.c_2.k562 <- processflow(exp.parp2_2.c_2.k562)

val.effectsizes <- data.frame(intgroup = c("BLM:FANCM:SPIDR", "BLM:FANCM:SPIDR",
                                           "PARP", "PARP",
                                           "PARP", "PARP",
                                           "MND1:PSMC3IP:RAD54L",
                                           "MND1:PSMC3IP:RAD54L",
                                           "MND1:PSMC3IP:RAD54L"),
                              int = c("BLM:SPIDR", "BLM:SPIDR", 
                                      "PARP1:RNASEH2C", "PARP1:RNASEH2C", 
                                      "PARP2:RNASEH2C", "PARP2:RNASEH2C", 
                                      "MND1:PSMC3IP", 
                                      "MND1:RAD54L", 
                                      "PSMC3IP:RAD54L"),
                              rep = c(rep(c("R1", "R2"), 3), "R1", "R1", "R1"),
                              D = c(flow.blm_1.sp_2.k562[[3]], flow.blm_2.sp_1.k562[[3]],
                                    flow.parp1_1.c_1.k562[[3]], flow.parp1_2.c_2.k562[[3]],
                                    flow.parp2_1.c_1.k562[[3]], flow.parp2_2.c_2.k562[[3]],
                                    flow.mnd_2.hop2_2.k562[[3]],
                                    flow.mnd_2.rad_2.k562[[3]],
                                    flow.hop2_2.rad_2.k562[[3]]))

val.mean.effect <- val.effectsizes %>% group_by(int) %>% summarise(D.mean = mean(D), diff = (max(D) - min(D))/2)

all.val.mean.all <- ints2 %>% group_by(int) %>%
  summarise(GI.quad.mean = mean(GI.z), GI.quad.sd = sd(GI.z),
            GI.add.mean = mean(GI.add.z), GI.add.sd = sd(GI.add.z))

val.screen.join <- inner_join(all.val.mean.all, val.mean.effect)

val.screen <- ggplot(val.screen.join, aes(GI.quad.mean, D.mean)) + 
  geom_vline(xintercept = 0,alpha = .3) +
  geom_hline(yintercept = 0, alpha = .3) +
  geom_abline(alpha = .3) +
  geom_point(size = 3) + 
  geom_text_repel(aes(label = int), size = 3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "none") +
  ylab("Validation Interaction Score") + xlab("Screen Interaction Score") +
  coord_fixed() +
  annotate(geom = "text", x = -10, y = 8, hjust = 0, 
           label = paste("r =",round(cor(val.screen.join$GI.quad.mean, val.screen.join$D.mean), 3)))

ggsave("../FIGURES/FIGURE-S3/validation_screen_intscore_compare.svg", val.screen,
       device = "svg", height = 4, width = 4)
