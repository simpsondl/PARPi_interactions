library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(outliers)

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
                  group = interaction(population, replicate),
                  linetype = replicate)) +
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
  
  return(list(flow, warnings, ttest, plt))
}


h2c <- file.path("../DATA/ValidationExperiments/rnaseh2c-nt-bfp_parp1-gfp.txt")
parp2 <- file.path("../DATA/ValidationExperiments/nt-parp1-bfp_rnaseh2c-gfp.txt")
triple <- file.path("../DATA/ValidationExperiments/rnaseh2c-parp2-bfp_parp1-gfp.txt")

h2c.res <- processflow(h2c)
parp2.res <- processflow(parp2)
triple.res <- processflow(triple)

all.res <- rbind(h2c.res[[1]], parp2.res[[1]], triple.res[[1]])
all.res$Exp <- rep(c("RNASEH2C", "PARP2", "TRIPLE"), each = 15)

parp1.means <- colMeans(all.res[all.res$population == "GFP+",3:(ncol(all.res)-1)])
sf <- all.res %>% group_by(Exp) %>% filter(population == "GFP+") %>% summarise(day6 = mean(day6),
                                                                               day9 = mean(day9),
                                                                               day12 = mean(day12),
                                                                               day15 = mean(day15),)
sf2 <- sf
sf2[,3:ncol(sf2)] <- sweep(sf2[,3:ncol(sf2)], 2, parp1.means[2:length(parp1.means)], "/")
sf2[,3:ncol(sf2)] <- 1/sf2[,3:ncol(sf2)]

all.res$day9.Adj[all.res$Exp == "RNASEH2C"] <- all.res$day9[all.res$Exp == "RNASEH2C"] * sf2$day9[sf2$Exp == "RNASEH2C"]
all.res$day9.Adj[all.res$Exp == "PARP2"] <- all.res$day9[all.res$Exp == "PARP2"] * sf2$day9[sf2$Exp == "PARP2"]
all.res$day9.Adj[all.res$Exp == "TRIPLE"] <- all.res$day9[all.res$Exp == "TRIPLE"] * sf2$day9[sf2$Exp == "TRIPLE"]

all.res$day12.Adj[all.res$Exp == "RNASEH2C"] <- all.res$day12[all.res$Exp == "RNASEH2C"] * sf2$day12[sf2$Exp == "RNASEH2C"]
all.res$day12.Adj[all.res$Exp == "PARP2"] <- all.res$day12[all.res$Exp == "PARP2"] * sf2$day12[sf2$Exp == "PARP2"]
all.res$day12.Adj[all.res$Exp == "TRIPLE"] <- all.res$day12[all.res$Exp == "TRIPLE"] * sf2$day12[sf2$Exp == "TRIPLE"]

all.res$day15.Adj[all.res$Exp == "RNASEH2C"] <- all.res$day15[all.res$Exp == "RNASEH2C"] * sf2$day15[sf2$Exp == "RNASEH2C"]
all.res$day15.Adj[all.res$Exp == "PARP2"] <- all.res$day15[all.res$Exp == "PARP2"] * sf2$day15[sf2$Exp == "PARP2"]
all.res$day15.Adj[all.res$Exp == "TRIPLE"] <- all.res$day15[all.res$Exp == "TRIPLE"] * sf2$day15[sf2$Exp == "TRIPLE"]


all.m <- melt(all.res[,c(1:2,7,3,9,10,8)], id.vars = c("population", "replicate", "Exp"))
all.m$day <- gsub("day", "", all.m$variable)
all.m$day <- gsub(".Adj", "", all.m$day)
all.m$day <- as.numeric(all.m$day)
# Remove the chaff
all.m2 <- all.m[all.m$population != "-/-" & all.m$population != "Exp",]
all.m2 <- all.m2[!(all.m2$population == "BFP+" & all.m2$Exp != "TRIPLE"),]

all.m2$population <- as.character(all.m2$population)
all.m2$population[all.m2$population == "GFP+"] <- "PARP1"
all.m2$population[all.m2$population == "BFP+" & all.m2$Exp == "RNASEH2C"] <- "RNASEH2C"
all.m2$population[all.m2$population == "BFP+" & all.m2$Exp == "PARP2"] <- "PARP2"
all.m2$population[all.m2$population == "BFP+GFP+" & all.m2$Exp == "RNASEH2C"] <- "PARP1:RNASEH2C"
all.m2$population[all.m2$population == "BFP+GFP+" & all.m2$Exp == "PARP2"] <- "PARP1:PARP2"
all.m2$population[all.m2$population == "BFP+" & all.m2$Exp == "TRIPLE"] <- "PARP2:RNASEH2C"
all.m2$population[all.m2$population == "BFP+GFP+" & all.m2$Exp == "TRIPLE"] <- "PARP1:PARP2:RNASEH2C"

all.m2$group <- paste(all.m2$population, all.m2$Exp, sep = ":")

rnaseh2c.diffs <- all.res[!(all.res$population %in% c("-/-", "Exp")) & all.res$Exp == "RNASEH2C",]
rnaseh2c.diffs2 <- rnaseh2c.diffs %>% group_by(replicate) %>% 
  summarise(Obs.day9 = day9.Adj[population == "BFP+GFP+"],
            Exp.day9 = day9.Adj[population == "BFP+"]+day9.Adj[population == "GFP+"],
            Obs.day12 = day12.Adj[population == "BFP+GFP+"],
            Exp.day12 = day12.Adj[population == "BFP+"]+day12.Adj[population == "GFP+"],
            Obs.day15 = day15.Adj[population == "BFP+GFP+"],
            Exp.day15 = day15.Adj[population == "BFP+"]+day15.Adj[population == "GFP+"])
rnaseh2c.diffs2$Eff.day9 <- rnaseh2c.diffs2$Obs.day9 - rnaseh2c.diffs2$Exp.day9
rnaseh2c.diffs2$Eff.day12 <- rnaseh2c.diffs2$Obs.day12 - rnaseh2c.diffs2$Exp.day12
rnaseh2c.diffs2$Eff.day15 <- rnaseh2c.diffs2$Obs.day15 - rnaseh2c.diffs2$Exp.day15

parp2.diffs <- all.res[!(all.res$population %in% c("-/-", "Exp")) & all.res$Exp == "PARP2",]
parp2.diffs2 <- parp2.diffs %>% group_by(replicate) %>% 
  summarise(Obs.day9 = day9.Adj[population == "BFP+GFP+"],
            Exp.day9 = day9.Adj[population == "BFP+"]+day9.Adj[population == "GFP+"],
            Obs.day12 = day12.Adj[population == "BFP+GFP+"],
            Exp.day12 = day12.Adj[population == "BFP+"]+day12.Adj[population == "GFP+"],
            Obs.day15 = day15.Adj[population == "BFP+GFP+"],
            Exp.day15 = day15.Adj[population == "BFP+"]+day15.Adj[population == "GFP+"])
parp2.diffs2$Eff.day9 <- parp2.diffs2$Obs.day9 - parp2.diffs2$Exp.day9
parp2.diffs2$Eff.day12 <- parp2.diffs2$Obs.day12 - parp2.diffs2$Exp.day12
parp2.diffs2$Eff.day15 <- parp2.diffs2$Obs.day15 - parp2.diffs2$Exp.day15

exp.trip.df <- data.frame(replicate = c("R1", "R2", "R3"),
                          day6 = 0,
                          day9 = c(min(all.res$day9.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     min(all.res$day9.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     min(rnaseh2c.diffs2$Eff.day9) + min(parp2.diffs2$Eff.day9),
                                   median(all.res$day9.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     median(all.res$day9.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     median(rnaseh2c.diffs2$Eff.day9) + median(parp2.diffs2$Eff.day9),
                                   max(all.res$day9.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     max(all.res$day9.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     max(rnaseh2c.diffs2$Eff.day9) + max(parp2.diffs2$Eff.day9)),
                          day12 = c(min(all.res$day12.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     min(all.res$day12.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     min(rnaseh2c.diffs2$Eff.day12) + min(parp2.diffs2$Eff.day12),
                                   median(all.res$day12.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     median(all.res$day12.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     median(rnaseh2c.diffs2$Eff.day12) + median(parp2.diffs2$Eff.day12),
                                   max(all.res$day12.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     max(all.res$day12.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     max(rnaseh2c.diffs2$Eff.day12) + max(parp2.diffs2$Eff.day12)),
                          day15 = c(min(all.res$day15.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     min(all.res$day15.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     min(rnaseh2c.diffs2$Eff.day15) + min(parp2.diffs2$Eff.day15),
                                   median(all.res$day15.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     median(all.res$day15.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     median(rnaseh2c.diffs2$Eff.day15) + median(parp2.diffs2$Eff.day15),
                                   max(all.res$day15.Adj[all.res$population == "GFP+" & all.res$Exp == "TRIPLE"]) + #PARP1 alone
                                     max(all.res$day15.Adj[all.res$population == "BFP+" & all.res$Exp == "TRIPLE"]) + #PARP2xRNASEH2C
                                     max(rnaseh2c.diffs2$Eff.day15) + max(parp2.diffs2$Eff.day15)))


exp.trip.m <- melt(exp.trip.df)
exp.trip.m$day <- gsub("day", "", exp.trip.m$variable)
exp.trip.m$day <- as.numeric(exp.trip.m$day)
exp.trip.m$population <- "Expected PARP1:PARP2:RNASEH2C"
exp.trip.m$Exp <- "Expected"
exp.trip.m$group <- "Expected"
exp.trip.m <- exp.trip.m[,colnames(all.m2)]

all.m2 <- rbind(all.m2, exp.trip.m)

all.m2$group <- factor(all.m2$group, levels = c("PARP2:PARP2", "RNASEH2C:RNASEH2C", 
                                                "PARP1:PARP2", "PARP1:RNASEH2C", "PARP1:TRIPLE",
                                                "PARP1:RNASEH2C:RNASEH2C", "PARP1:PARP2:PARP2",
                                                "PARP2:RNASEH2C:TRIPLE", "Expected",
                                                "PARP1:PARP2:RNASEH2C:TRIPLE"))
all.m2$population <- factor(all.m2$population, levels = c("PARP1", "PARP2", "RNASEH2C", 
                                                          "PARP1:RNASEH2C", "PARP1:PARP2", 
                                                          "PARP2:RNASEH2C", "Expected PARP1:PARP2:RNASEH2C",
                                                          "PARP1:PARP2:RNASEH2C"))

errors <- all.res %>% group_by(population, Exp) %>% 
  summarise(day9.mean = mean(day9.Adj), day9.sd = sd(day9.Adj),
            day12.mean = mean(day12.Adj), day12.sd = sd(day12.Adj),
            day15.mean = mean(day15.Adj), day15.sd = sd(day15.Adj))

errors.df <- data.frame(population = errors$population ,
                        Exp = errors$Exp,
                        day9min = errors$day9.mean - errors$day9.sd,
                        day9mean = errors$day9.mean,
                        day9max = errors$day9.mean + errors$day9.sd,
                        day12min = errors$day12.mean - errors$day12.sd,
                        day12mean = errors$day12.mean,
                        day12max = errors$day12.mean + errors$day12.sd,
                        day15min = errors$day15.mean - errors$day15.sd,
                        day15mean = errors$day15.mean,
                        day15max = errors$day15.mean + errors$day15.sd)

errors.df.m <- melt(errors.df, id.vars = c("population", "Exp"))
errors.df.m$day <- gsub("day", "", errors.df.m$variable)
errors.df.m$stat <- gsub("[0-9]+", "", errors.df.m$day)
errors.df.m$day <- as.numeric(gsub("[a-z]+", "", errors.df.m$day))

errors.wide <- reshape(data = errors.df.m[,c(1,2,4:6)], idvar = c("population", "Exp", "day"), 
                       timevar = "stat", direction = "wide")
# Remove the chaff
errors.wide <- errors.wide[errors.wide$population != "-/-" & errors.wide$population != "Exp",]
errors.wide <- errors.wide[!(errors.wide$population == "BFP+" & errors.wide$Exp != "TRIPLE"),]

errors.wide$population <- as.character(errors.wide$population)
errors.wide$population[errors.wide$population == "GFP+"] <- "PARP1"
errors.wide$population[errors.wide$population == "BFP+GFP+" & errors.wide$Exp == "RNASEH2C"] <- "PARP1:RNASEH2C"
errors.wide$population[errors.wide$population == "BFP+GFP+" & errors.wide$Exp == "PARP2"] <- "PARP1:PARP2"
errors.wide$population[errors.wide$population == "BFP+" & errors.wide$Exp == "TRIPLE"] <- "PARP2:RNASEH2C"
errors.wide$population[errors.wide$population == "BFP+GFP+" & errors.wide$Exp == "TRIPLE"] <- "PARP1:PARP2:RNASEH2C"

ggplot(all.m2) + 
  geom_point(aes(day, value, col = group)) + 
  #geom_errorbar(data = errors.wide, aes(x = day, ymin = value.min, ymax = value.max)) +
  stat_summary(aes(day, value, col = group), geom = "line", fun = mean, size = 2) +
  scale_color_manual(values = c("#117733", "#CC6677", "#ddcc77","#ddcc77","#ddcc77","#bb2cbb", "#861d86","#500d50", "#999999", "black"),
                     #labels = c("PARP1", "PARP1:PARP2", "PARP1", "PARP1:RNASEH2C", "PARP1", "PARP1:PARP2:RNASEH2C", "PARP2:RNASEH2C"),
                     name = "Interaction") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"))
