library(readr)
library(dplyr)
library(ggplot2)

# Load functions
source("helper_functions.R")

con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_avg.txt", sep = ""))
con.tau.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_tau_oi_avg.txt", sep = ""))
con.id.map <- read_tsv("../DATA/pseudogenecombination_id_map.txt")
gene.id.map <- read_tsv("../DATA/genecombination_id_map.txt")

con.oi.avg <- compute_construct_differential_scores(con.gamma.oi.avg, con.tau.oi.avg)

gamma.gene.gi <- compute_gene_interaction_scores(con.gamma.oi.avg, "GI.z")
tau.gene.gi <- compute_gene_interaction_scores(con.tau.oi.avg, "GI.z")
nu.gene.gi <- compute_gene_interaction_scores(con.oi.avg, "DGI.z")

gamma.bound.sig <- 1.524 #qvalue <= .05
tau.bound.sig <- 1.527 #qvalue <= .05
nu.bound.sig <- 2.975 #qvalue <= .20

gamma.bound.hc <- 2.723 #local FDR <= .01
tau.bound.hc <- 2.853 #local FDR <= .01
nu.bound.hc <- 2.975 #qvalue <= .20

gi.cmp <-inner_join(gamma.gene.gi[,1:3], tau.gene.gi[,1:3], by = c("GeneCombinationID", "Category"), suffix = c(".Gamma", ".Tau"))
gi.cmp <-inner_join(gi.cmp, nu.gene.gi[,1:3], by = c("GeneCombinationID", "Category"))
colnames(gi.cmp)[3:5] <- c("Gamma", "Tau", "Nu")

gi.cmp2 <- inner_join(gi.cmp, gene.id.map, 
                      by = c("GeneCombinationID")) %>%
  filter(Category == "X+Y")

gi.cmp2$InteractionCategory <- NA
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & abs(gi.cmp2$Tau) < tau.bound.sig] <- "No interaction"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma > 0 & gi.cmp2$Tau > 0] <- "Interaction, unmodified - buffering"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma < 0 & gi.cmp2$Tau < 0] <- "Interaction, unmodified - synthetic lethal"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig & gi.cmp2$Gamma < 0 & gi.cmp2$Tau > 0] <- "Interaction, reversed - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig & gi.cmp2$Gamma > 0 & gi.cmp2$Tau < 0] <- "Interaction, reversed - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma < 0 & gi.cmp2$Tau > 0] <- "Interaction, weakly reversed - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) >= gamma.bound.sig & abs(gi.cmp2$Tau) >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig & gi.cmp2$Gamma > 0 & gi.cmp2$Tau < 0] <- "Interaction, weakly reversed - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau >= tau.bound.sig & gi.cmp2$Nu >= nu.bound.sig] <- "Interaction, novel - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau <= -1*tau.bound.sig & gi.cmp2$Nu <= -1*nu.bound.sig] <- "Interaction, novel - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau >= tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, weakly novel - nu positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Gamma) < gamma.bound.sig & gi.cmp2$Tau <= -1*tau.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, weakly novel - nu negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig] <- "Interaction, masked - gamma positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & abs(gi.cmp2$Nu) >= nu.bound.sig] <- "Interaction, masked - gamma negative"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, weakly masked - gamma positive"
gi.cmp2$InteractionCategory[abs(gi.cmp2$Tau) < tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & abs(gi.cmp2$Nu) < nu.bound.sig] <- "Interaction, weakly masked - gamma negative"
gi.cmp2$InteractionCategory[gi.cmp2$Tau >= tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & gi.cmp2$Nu >= nu.bound.sig] <- "Interaction, modified - buffering, nu positive"
gi.cmp2$InteractionCategory[gi.cmp2$Tau >= tau.bound.sig & gi.cmp2$Gamma >= gamma.bound.sig & gi.cmp2$Nu <= -1*nu.bound.sig] <- "Interaction, modified - buffering, nu negative"
gi.cmp2$InteractionCategory[gi.cmp2$Tau <= -1*tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & gi.cmp2$Nu >= nu.bound.sig] <- "Interaction, modified - synthetic lethal, nu positive"
gi.cmp2$InteractionCategory[gi.cmp2$Tau <= -1*tau.bound.sig & gi.cmp2$Gamma <= -1*gamma.bound.sig & gi.cmp2$Nu <= -1*nu.bound.sig] <- "Interaction, modified - synthetic lethal, nu negative"

write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "No interaction",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_nointeraction.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, unmodified - buffering",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_unmodified_buffering.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, unmodified - synthetic lethal",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_unmodified_synlethal.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, reversed - nu positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_reversed_nupositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, reversed - nu negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_reversed_nunegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, weakly reversed - nu positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_weaklyreversed_nupositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, weakly reversed - nu negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_weaklyreversed_nunegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, novel - nu positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_novel_nupositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, novel - nu negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_novel_nunegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, weakly novel - nu positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_weaklynovel_nupositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, weakly novel - nu negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_weaklynovel_nunegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, masked - gamma positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_masked_gammapositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, masked - gamma negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_masked_gammanegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, weakly masked - gamma positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_weaklymasked_gammapositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, weakly masked - gamma negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_weaklymasked_gammanegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, modified - buffering, nu positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_modified_buffering_nupositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, modified - buffering, nu negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_modified_buffering_nunegative.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, modified - synthetic lethal, nu positive",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_modified_synlethal_nupositive.txt")
write_tsv(gi.cmp2[gi.cmp2$InteractionCategory == "Interaction, modified - synthetic lethal, nu negative",c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_modified_synlethal_nunegative.txt")
write_tsv(gi.cmp2[,c(6:8,1:5,9)], 
          "../DATA/Interaction_Scores/InteractionCategories/genecombination_interactioncategories_all_interactions.txt")



counts <- data.frame(Gene = unique(c(gi.cmp2$Gene2, gi.cmp2$Gene1)), 
                     M = sum(gi.cmp2$Category == "X+Y"),
                     m.unmodified = sum(grepl("Interaction, unmodified", gi.cmp2$InteractionCategory)),
                     m.unmodbuff = sum(gi.cmp2$InteractionCategory == "Interaction, unmodified - buffering"),
                     m.unmodsynlethal = sum(gi.cmp2$InteractionCategory == "Interaction, unmodified - synthetic lethal"),
                     m.modified = sum(grepl("Interaction, modified", gi.cmp2$InteractionCategory)),
                     m.modbuffpos = sum(gi.cmp2$InteractionCategory == "Interaction, modified - buffering, nu positive"),
                     m.modbuffneg = sum(gi.cmp2$InteractionCategory == "Interaction, modified - buffering, nu negative"),
                     m.modsynlethpos = sum(gi.cmp2$InteractionCategory == "Interaction, modified - synthetic lethal, nu positive"),
                     m.modsynlethneg = sum(gi.cmp2$InteractionCategory == "Interaction, modified - synthetic lethal, nu negative"),
                     m.masked = sum(grepl("Interaction, masked", gi.cmp2$InteractionCategory)),
                     m.maskedpos = sum(gi.cmp2$InteractionCategory == "Interaction, masked - gamma positive"),
                     m.maskedneg = sum(gi.cmp2$InteractionCategory == "Interaction, masked - gamma negative"),
                     m.weakmasked = sum(grepl("Interaction, weakly masked", gi.cmp2$InteractionCategory)),
                     m.weakmaskedpos = sum(gi.cmp2$InteractionCategory == "Interaction, weakly masked - gamma positive"),
                     m.weakmaskedneg = sum(gi.cmp2$InteractionCategory == "Interaction, weakly masked - gamma negative"),
                     m.novel = sum(grepl("Interaction, novel", gi.cmp2$InteractionCategory)),
                     m.novelpos = sum(gi.cmp2$InteractionCategory == "Interaction, novel - nu positive"),
                     m.novelneg = sum(gi.cmp2$InteractionCategory == "Interaction, novel - nu negative"),
                     m.weaknovel = sum(grepl("Interaction, weakly novel", gi.cmp2$InteractionCategory)),
                     m.weaknovelpos = sum(gi.cmp2$InteractionCategory == "Interaction, weakly novel - nu positive"),
                     m.weaknovelneg = sum(gi.cmp2$InteractionCategory == "Interaction, weakly novel - nu negative"),
                     m.reversed = sum(grepl("Interaction, reversed", gi.cmp2$InteractionCategory)),
                     m.reversedpos = sum(gi.cmp2$InteractionCategory == "Interaction, reversed - nu positive"),
                     m.reversedneg = sum(gi.cmp2$InteractionCategory == "Interaction, reversed - nu negative"),
                     m.weakreversed = sum(grepl("Interaction, weakly reversed", gi.cmp2$InteractionCategory)),
                     m.weakreversedpos = sum(gi.cmp2$InteractionCategory == "Interaction, weakly reversed - nu positive"),
                     m.weakreversedneg = sum(gi.cmp2$InteractionCategory == "Interaction, weakly reversed - nu negative"))

counts$N <- NA
counts$n.unmodified <- NA
counts$n.unmodbuff <- NA
counts$n.unmodsynlethal <- NA
counts$n.modified <- NA
counts$n.modbuffpos <- NA
counts$n.modbuffneg <- NA 
counts$n.modsynlethpos <- NA
counts$n.modsynlethneg <- NA
counts$n.masked <- NA
counts$n.maskedpos<- NA
counts$n.maskedneg<- NA
counts$n.weakmasked<- NA
counts$n.weakmaskedpos <- NA
counts$n.weakmaskedneg<- NA
counts$n.novel<- NA
counts$n.novelpos<- NA
counts$n.novelneg <- NA
counts$n.weaknovel <- NA
counts$n.weaknovelpos <- NA
counts$n.weaknovelneg <- NA
counts$n.reversed <- NA
counts$n.reversedpos<- NA
counts$n.reversedneg <- NA
counts$n.weakreversed<- NA
counts$n.weakreversedpos<- NA
counts$n.weakreversedneg<- NA


for(i in counts$Gene){
  counts$N[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName))
  counts$n.unmodified[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, unmodified", gi.cmp2$InteractionCategory))
  counts$n.unmodbuff[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, unmodified - buffering")
  counts$n.unmodsynlethal[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, unmodified - synthetic lethal")
  counts$n.modified[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, modified", gi.cmp2$InteractionCategory))
  counts$n.modbuffpos[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, modified - buffering, nu positive")
  counts$n.modbuffneg[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, modified - buffering, nu negative") 
  counts$n.modsynlethpos[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, modified - synthetic lethal, nu positive")
  counts$n.modsynlethneg[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, modified - synthetic lethal, nu negative")
  counts$n.masked[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, masked", gi.cmp2$InteractionCategory))
  counts$n.maskedpos[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, masked - gamma positive")
  counts$n.maskedneg[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, masked - gamma negative")
  counts$n.weakmasked[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, weakly masked", gi.cmp2$InteractionCategory))
  counts$n.weakmaskedpos[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, weakly masked - gamma positive")
  counts$n.weakmaskedneg[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, weakly masked - gamma negative")
  counts$n.novel[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, novel", gi.cmp2$InteractionCategory))
  counts$n.novelpos[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, novel - nu positive")
  counts$n.novelneg[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, novel - nu negative")
  counts$n.weaknovel[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, weakly novel", gi.cmp2$InteractionCategory))
  counts$n.weaknovelpos[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, weakly novel - nu positive")
  counts$n.weaknovelneg[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, weakly novel - nu negative")
  counts$n.reversed[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, reversed", gi.cmp2$InteractionCategory))
  counts$n.reversedpos[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, reversed - nu positive")
  counts$n.reversedneg[counts$Gene == i] <- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, reversed - nu negative")
  counts$n.weakreversed[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & grepl("Interaction, weakly reversed", gi.cmp2$InteractionCategory))
  counts$n.weakreversedpos[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, weakly reversed - nu positive")
  counts$n.weakreversedneg[counts$Gene == i]<- sum(grepl(paste0("^", i, ":|:", i, "$"),gi.cmp2$GeneCombinationName) & gi.cmp2$InteractionCategory == "Interaction, weakly reversed - nu negative")
}

write_tsv(counts, "../DATA/Interaction_Scores/InteractionCategories/gene_interactioncategory_counts.txt")

res <- data.frame(Gene = unique(c(gi.cmp2$Gene2, gi.cmp2$Gene1)),
                     Unmodified = NA,
                     Unmodbuff = NA,
                     Unmodsynlethal = NA,
                     Modified = NA,
                     Modbuffpos = NA,
                     Modbuffneg = NA,
                     Modsynlethpos = NA,
                     Modsynlethneg = NA,
                     Masked = NA,
                     Maskedpos = NA,
                     Maskedneg = NA,
                     Weakmasked = NA,
                     Weakmaskedpos = NA,
                     Weakmaskedneg = NA,
                     Novel = NA,
                     Novelpos = NA,
                     Novelneg = NA,
                     Weaknovel = NA,
                     Weaknovelpos = NA,
                     Weaknovelneg = NA,
                     Reversed = NA,
                     Reversedpos = NA,
                     Reversedneg = NA,
                     Weakreversed = NA,
                     Weakreversedpos = NA,
                     Weakreversedneg = NA)


for( i in res$Gene){
  
  res$Unmodified[res$Gene == i] <- fisher.test(matrix(c(counts$n.unmodified[counts$Gene == i], 
                                                        counts$N[counts$Gene == i] - counts$n.unmodified[counts$Gene == i], 
                                                        counts$m.unmodified[counts$Gene == i] - counts$n.unmodified[counts$Gene == i], 
                                                        counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                          counts$m.unmodified[counts$Gene == i] + counts$n.unmodified[counts$Gene == i]), 2))$p.value
    res$Unmodbuff[res$Gene == i] <- fisher.test(matrix(c(counts$n.unmodbuff[counts$Gene == i], 
                                                         counts$N[counts$Gene == i] - counts$n.unmodbuff[counts$Gene == i], 
                                                         counts$m.unmodbuff[counts$Gene == i] - counts$n.unmodbuff[counts$Gene == i], 
                                                         counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                           counts$m.unmodbuff[counts$Gene == i] + counts$n.unmodbuff[counts$Gene == i]), 2))$p.value
    res$Unmodsynlethal[res$Gene == i] <- fisher.test(matrix(c(counts$n.unmodsynlethal[counts$Gene == i], 
                                                              counts$N[counts$Gene == i] - counts$n.unmodsynlethal[counts$Gene == i], 
                                                              counts$m.unmodsynlethal[counts$Gene == i] - counts$n.unmodsynlethal[counts$Gene == i], 
                                                              counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                                counts$m.unmodsynlethal[counts$Gene == i] + counts$n.unmodsynlethal[counts$Gene == i]), 2))$p.value
    res$Modified[res$Gene == i] <- fisher.test(matrix(c(counts$n.modified[counts$Gene == i], 
                                                        counts$N[counts$Gene == i] - counts$n.modified[counts$Gene == i], 
                                                        counts$m.modified[counts$Gene == i] - counts$n.modified[counts$Gene == i], 
                                                        counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                          counts$m.modified[counts$Gene == i] + counts$n.modified[counts$Gene == i]), 2))$p.value
    res$Modbuffpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.modbuffpos[counts$Gene == i], 
                                                          counts$N[counts$Gene == i] - counts$n.modbuffpos[counts$Gene == i], 
                                                          counts$m.modbuffpos[counts$Gene == i] - counts$n.modbuffpos[counts$Gene == i], 
                                                          counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                            counts$m.modbuffpos[counts$Gene == i] + counts$n.modbuffpos[counts$Gene == i]), 2))$p.value
    res$Modbuffneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.modbuffneg[counts$Gene == i], 
                                                          counts$N[counts$Gene == i] - counts$n.modbuffneg[counts$Gene == i], 
                                                          counts$m.modbuffneg[counts$Gene == i] - counts$n.modbuffneg[counts$Gene == i], 
                                                          counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                            counts$m.modbuffneg[counts$Gene == i] + counts$n.modbuffneg[counts$Gene == i]), 2))$p.value
    res$Modsynlethpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.modsynlethpos[counts$Gene == i], 
                                                             counts$N[counts$Gene == i] - counts$n.modsynlethpos[counts$Gene == i], 
                                                             counts$m.modsynlethpos[counts$Gene == i] - counts$n.modsynlethpos[counts$Gene == i], 
                                                             counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                               counts$m.modsynlethpos[counts$Gene == i] + counts$n.modsynlethpos[counts$Gene == i]), 2))$p.value
    res$Modsynlethneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.modsynlethneg[counts$Gene == i], 
                                                             counts$N[counts$Gene == i] - counts$n.modsynlethneg[counts$Gene == i], 
                                                             counts$m.modsynlethneg[counts$Gene == i] - counts$n.modsynlethneg[counts$Gene == i], 
                                                             counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                               counts$m.modsynlethneg[counts$Gene == i] + counts$n.modsynlethneg[counts$Gene == i]), 2))$p.value
    res$Masked[res$Gene == i] <- fisher.test(matrix(c(counts$n.masked[counts$Gene == i], 
                                                      counts$N[counts$Gene == i] - counts$n.masked[counts$Gene == i], 
                                                      counts$m.masked[counts$Gene == i] - counts$n.masked[counts$Gene == i], 
                                                      counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                        counts$m.masked[counts$Gene == i] + counts$n.masked[counts$Gene == i]), 2))$p.value
    res$Maskedpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.maskedpos[counts$Gene == i], 
                                                         counts$N[counts$Gene == i] - counts$n.maskedpos[counts$Gene == i], 
                                                         counts$m.maskedpos[counts$Gene == i] - counts$n.maskedpos[counts$Gene == i], 
                                                         counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                           counts$m.maskedpos[counts$Gene == i] + counts$n.maskedpos[counts$Gene == i]), 2))$p.value
    res$Maskedneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.maskedneg[counts$Gene == i], 
                                                         counts$N[counts$Gene == i] - counts$n.maskedneg[counts$Gene == i], 
                                                         counts$m.maskedneg[counts$Gene == i] - counts$n.maskedneg[counts$Gene == i], 
                                                         counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                           counts$m.maskedneg[counts$Gene == i] + counts$n.maskedneg[counts$Gene == i]), 2))$p.value
    res$Weakmasked[res$Gene == i] <- fisher.test(matrix(c(counts$n.weakmasked[counts$Gene == i], 
                                                          counts$N[counts$Gene == i] - counts$n.weakmasked[counts$Gene == i], 
                                                          counts$m.weakmasked[counts$Gene == i] - counts$n.weakmasked[counts$Gene == i], 
                                                          counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                            counts$m.weakmasked[counts$Gene == i] + counts$n.weakmasked[counts$Gene == i]), 2))$p.value
    res$Weakmaskedpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.weakmaskedpos[counts$Gene == i], 
                                                             counts$N[counts$Gene == i] - counts$n.weakmaskedpos[counts$Gene == i], 
                                                             counts$m.weakmaskedpos[counts$Gene == i] - counts$n.weakmaskedpos[counts$Gene == i], 
                                                             counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                               counts$m.weakmaskedpos[counts$Gene == i] + counts$n.weakmaskedpos[counts$Gene == i]), 2))$p.value
    res$Weakmaskedneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.weakmaskedneg[counts$Gene == i], 
                                                             counts$N[counts$Gene == i] - counts$n.weakmaskedneg[counts$Gene == i], 
                                                             counts$m.weakmaskedneg[counts$Gene == i] - counts$n.weakmaskedneg[counts$Gene == i], 
                                                             counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                               counts$m.weakmaskedneg[counts$Gene == i] + counts$n.weakmaskedneg[counts$Gene == i]), 2))$p.value
    res$Novel[res$Gene == i] <- fisher.test(matrix(c(counts$n.novel[counts$Gene == i], 
                                                     counts$N[counts$Gene == i] - counts$n.novel[counts$Gene == i], 
                                                     counts$m.novel[counts$Gene == i] - counts$n.novel[counts$Gene == i], 
                                                     counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                       counts$m.novel[counts$Gene == i] + counts$n.novel[counts$Gene == i]), 2))$p.value
    res$Novelpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.novelpos[counts$Gene == i], 
                                                        counts$N[counts$Gene == i] - counts$n.novelpos[counts$Gene == i], 
                                                        counts$m.novelpos[counts$Gene == i] - counts$n.novelpos[counts$Gene == i], 
                                                        counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                          counts$m.novelpos[counts$Gene == i] + counts$n.novelpos[counts$Gene == i]), 2))$p.value
    res$Novelneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.novelneg[counts$Gene == i], 
                                                        counts$N[counts$Gene == i] - counts$n.novelneg[counts$Gene == i], 
                                                        counts$m.novelneg[counts$Gene == i] - counts$n.novelneg[counts$Gene == i], 
                                                        counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                          counts$m.novelneg[counts$Gene == i] + counts$n.novelneg[counts$Gene == i]), 2))$p.value
    res$Weaknovel[res$Gene == i] <- fisher.test(matrix(c(counts$n.weaknovel[counts$Gene == i], 
                                                         counts$N[counts$Gene == i] - counts$n.weaknovel[counts$Gene == i], 
                                                         counts$m.weaknovel[counts$Gene == i] - counts$n.weaknovel[counts$Gene == i], 
                                                         counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                           counts$m.weaknovel[counts$Gene == i] + counts$n.weaknovel[counts$Gene == i]), 2))$p.value
    res$Weaknovelpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.weaknovelpos[counts$Gene == i], 
                                                            counts$N[counts$Gene == i] - counts$n.weaknovelpos[counts$Gene == i], 
                                                            counts$m.weaknovelpos[counts$Gene == i] - counts$n.weaknovelpos[counts$Gene == i], 
                                                            counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                              counts$m.weaknovelpos[counts$Gene == i] + counts$n.weaknovelpos[counts$Gene == i]), 2))$p.value
    res$Weaknovelneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.weaknovelneg[counts$Gene == i], 
                                                            counts$N[counts$Gene == i] - counts$n.weaknovelneg[counts$Gene == i], 
                                                            counts$m.weaknovelneg[counts$Gene == i] - counts$n.weaknovelneg[counts$Gene == i], 
                                                            counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                              counts$m.weaknovelneg[counts$Gene == i] + counts$n.weaknovelneg[counts$Gene == i]), 2))$p.value
    res$Reversed[res$Gene == i] <- fisher.test(matrix(c(counts$n.reversed[counts$Gene == i], 
                                                        counts$N[counts$Gene == i] - counts$n.reversed[counts$Gene == i], 
                                                        counts$m.reversed[counts$Gene == i] - counts$n.reversed[counts$Gene == i], 
                                                        counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                          counts$m.reversed[counts$Gene == i] + counts$n.reversed[counts$Gene == i]), 2))$p.value
    res$Reversedpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.reversedpos[counts$Gene == i], 
                                                           counts$N[counts$Gene == i] - counts$n.reversedpos[counts$Gene == i], 
                                                           counts$m.reversedpos[counts$Gene == i] - counts$n.reversedpos[counts$Gene == i], 
                                                           counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                             counts$m.reversedpos[counts$Gene == i] + counts$n.reversedpos[counts$Gene == i]), 2))$p.value
    res$Reversedneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.reversedneg[counts$Gene == i], 
                                                           counts$N[counts$Gene == i] - counts$n.reversedneg[counts$Gene == i], 
                                                           counts$m.reversedneg[counts$Gene == i] - counts$n.reversedneg[counts$Gene == i], 
                                                           counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                             counts$m.reversedneg[counts$Gene == i] + counts$n.reversedneg[counts$Gene == i]), 2))$p.value
    res$Weakreversed[res$Gene == i] <- fisher.test(matrix(c(counts$n.weakreversed[counts$Gene == i], 
                                                            counts$N[counts$Gene == i] - counts$n.weakreversed[counts$Gene == i], 
                                                            counts$m.weakreversed[counts$Gene == i] - counts$n.weakreversed[counts$Gene == i], 
                                                            counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                              counts$m.weakreversed[counts$Gene == i] + counts$n.weakreversed[counts$Gene == i]), 2))$p.value
    res$Weakreversedpos[res$Gene == i] <- fisher.test(matrix(c(counts$n.weakreversedpos[counts$Gene == i], 
                                                               counts$N[counts$Gene == i] - counts$n.weakreversedpos[counts$Gene == i], 
                                                               counts$m.weakreversedpos[counts$Gene == i] - counts$n.weakreversedpos[counts$Gene == i], 
                                                               counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                                 counts$m.weakreversedpos[counts$Gene == i] + counts$n.weakreversedpos[counts$Gene == i]), 2))$p.value
    res$Weakreversedneg[res$Gene == i] <- fisher.test(matrix(c(counts$n.weakreversedneg[counts$Gene == i], 
                                                               counts$N[counts$Gene == i] - counts$n.weakreversedneg[counts$Gene == i], 
                                                               counts$m.weakreversedneg[counts$Gene == i] - counts$n.weakreversedneg[counts$Gene == i], 
                                                               counts$M[counts$Gene == i] - counts$N[counts$Gene == i] - 
                                                                 counts$m.weakreversedneg[counts$Gene == i] + counts$n.weakreversedneg[counts$Gene == i]), 2))$p.value
  
}

write_tsv(res, "../DATA/Interaction_Scores/InteractionCategories/gene_interactioncategory_enrichments.txt")

qunmodified <- qvalue(res$Unmodified)
qunmodbuff <- qvalue(res$Unmodbuff)
qunmodsynlethal <- qvalue(res$Unmodsynlethal)
qmodified <- qvalue(res$Modified)
qmodbuffpos <- qvalue(res$Modbuffpos)
#correcting floating point issue
res$Modbuffneg[res$Modbuffneg > .1] <- 1
qmodbuffneg <- qvalue(res$Modbuffneg)
qmodsynlethpos <- qvalue(res$Modsynlethpos)
qmodsynlethneg <- qvalue(res$Modsynlethneg)

qmasked <- qvalue(res$Masked)
qmaskedpos <- qvalue(res$Maskedpos)
qmaskedneg <- qvalue(res$Maskedneg)
qwmasked <- qvalue(res$Weakmasked)
qwmaskedpos <- qvalue(res$Weakmaskedpos)
qwmaskedneg <- qvalue(res$Weakmaskedneg)

qnovel <- qvalue(res$Novel)
qnovelpos <- qvalue(res$Novelpos)
qnovelneg <- qvalue(res$Novelneg)
qwnovel <- qvalue(res$Weaknovel)
qwnovelpos <- qvalue(res$Weaknovelpos)
qwnovelneg <- qvalue(res$Weaknovelneg)

qreverse <- qvalue(res$Reversed)
qreversepos <- qvalue(res$Reversedpos)
qreverseneg <- qvalue(res$Reversedneg)
# correct fp
res$Weakreversed[res$Weakreversed > .2] <- 1
qwreverse <- qvalue(res$Weakreversed)
qwreversepos <- qvalue(res$Weakreversedpos)
qwreverseneg <- qvalue(res$Weakreversedneg)

qs <- data.frame(Gene = res$Gene,
                 Q.Unmodified = qunmodified$qvalues,
                 Q.Unmodbuff = qunmodbuff$qvalues,
                 Q.Unmodsynlethal = qunmodsynlethal$qvalues,
                 Q.Modified = qmodified$qvalues,
                 Q.Modbuffpos = qmodbuffpos$qvalues,
                 Q.Modbuffneg = qmodbuffneg$qvalues,
                 Q.Modsynlethalpos = qmodsynlethpos$qvalues,
                 Q.Modsynlethalneg = qmodsynlethneg$qvalues,
                 Q.Masked = qmasked$qvalues,
                 Q.Maskedpos = qmaskedpos$qvalues,
                 Q.Maskedneg = qmaskedneg$qvalues,
                 Q.WeakMasked = qwmasked$qvalues,
                 Q.WeakMaskedpos = qwmaskedpos$qvalues,
                 Q.WeakMaskedneg = qwmaskedneg$qvalues,
                 Q.Novel = qnovel$qvalues,
                 Q.Novelpos = qnovelpos$qvalues,
                 Q.Novelneg = qnovelneg$qvalues,
                 Q.WeakNovel = qwnovel$qvalues,
                 Q.WeakNovelpos = qwnovelpos$qvalues,
                 Q.WeakNovelneg = qwnovelneg$qvalues,
                 Q.Reversed = qreverse$qvalues,
                 Q.Reversedpos = qreversepos$qvalues,
                 Q.Reversedneg = qreverseneg$qvalues,
                 Q.WeakReversed = qwreverse$qvalues,
                 Q.WeakReversedpos = qwreversepos$qvalues,
                 Q.WeakReversedneg = qwreverseneg$qvalues)

write_tsv(qs, "../DATA/Interaction_Scores/InteractionCategories/gene_interactioncategory_qvalues.txt")

gi.cmp2$Highlight <- NA
gi.cmp2$Highlight[gi.cmp2$InteractionScore.Nu >= nu.bound.sig] <- "Resistance Interaction"
gi.cmp2$Highlight[gi.cmp2$InteractionScore.Nu <= -1*nu.bound.sig] <- "Sensitizing Interaction"
gi.cmp2$Highlight[is.na(gi.cmp2$Highlight)] <- "No Interaction"
gi.cmp2$Highlight <- factor(gi.cmp2$Highlight, levels = c("No Interaction",
                                                          "Sensitizing Interaction", "Resistance Interaction"))

gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]
gi.samp <- gi.cmp2[sample(1:nrow(gi.cmp2),10000),]

nolabel_theme <-   theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.title.y = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.grid = element_blank(),
                         axis.ticks.length = unit(.2, "cm"),
                         legend.position = "none")
nolabel_axes <- list(xlab(""), ylab(""), ggtitle("")) 

plt <- ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau , col = Highlight)) +
        geom_hline(yintercept = 0, alpha = 0.5) +
        geom_hline(yintercept = gamma.bound.sig, alpha = 0.5, linetype = "dotted") +
        geom_hline(yintercept = -1*gamma.bound.sig, alpha = 0.5, linetype = "dotted") +
        geom_vline(xintercept = 0, alpha = 0.5) +
        geom_vline(xintercept = tau.bound.sig, alpha = 0.5, linetype = "dotted") +
        geom_vline(xintercept = -1*(tau.bound.sig), alpha= 0.5, linetype= "dotted") +
        geom_abline(slope = 1, alpha = 0.5) +
        geom_abline(slope = 1, intercept = nu.bound.sig, alpha = 0.5, linetype = "dotdash") +
        geom_abline(slope = 1, intercept = -1*nu.bound.sig, alpha = 0.5, linetype = "dotdash") +
        geom_point(size = .8, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "None") +
  xlab("Gamma IS") + ylab("Tau IS") +
  scale_x_continuous(breaks = c(-12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5,10)) +
  scale_y_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20)) +
        scale_color_manual(values = c("#999999", "#33716B", "darkorchid4" ))

plt.lab <- ggplot(gi.samp, aes(InteractionScore.Gamma, InteractionScore.Tau, col = Highlight)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = gamma.bound.sig, alpha = 0.5, linetype = "dotted") +
  geom_hline(yintercept = -1*gamma.bound.sig, alpha = 0.5, linetype = "dotted") +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = tau.bound.sig, alpha = 0.5, linetype = "dotted") +
  geom_vline(xintercept = -1*(tau.bound.sig), alpha= 0.5, linetype= "dotted") +
  geom_abline(slope = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 1, alpha = 0.5, linetype = "dotted") +
  geom_abline(slope = 1, intercept = -1, alpha = 0.5, linetype = "dotted") +
  geom_abline(slope = 1, intercept = nu.bound.sig, alpha = 0.5, linetype = "dotdash") +
  geom_abline(slope = 1, intercept = -1*nu.bound.sig, alpha = 0.5, linetype = "dotdash") +
  geom_abline(slope = 1, intercept = 2, alpha = 0.5, linetype = "dotted") +
  geom_abline(slope = 1, intercept = -2, alpha = 0.5, linetype = "dotted") +
  geom_point(size = .5, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Gamma IS") + ylab("Tau IS") +
  scale_x_continuous(breaks = c(-12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5,10)) +
  scale_y_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20)) +
  scale_color_manual(values = c("#999999", "#33716B", "#d8b365" )) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  labs(col = "Interaction Category")

ggsave("../FIGURES/FIGURE-4/gamma_tau_interaction_category_scatterplot.png", plt, device = "png",
       height = 8, width = 8, dpi = 300)
ggsave("../FIGURES/FIGURE-4/gamma_tau_interaction_category_scatterplot_labels.svg", plt.lab, device = "svg",
       height = 8, width = 10, dpi = 300)
