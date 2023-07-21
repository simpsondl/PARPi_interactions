library(readr)
library(dplyr)
library(ggplot2)

gtplot <- function(gene, gis){
  gis$Highlight <- gis$Gene1 == gene | gis$Gene2 == gene
  gis$Non <- NA
  gis$Non[gis$Highlight] <- gis$Gene1[gis$Highlight]
  gis$Non[!is.na(gis$Non) & gis$Non == gene] <- gis$Gene2[!is.na(gis$Non) & gis$Non == gene]
  gis <- left_join(gis, 
                   rho.single.gene.pheno,
                   by = c("Non" = "gene"))
  if(gene != "PARP1"){
    gis$rho[gis$Non == "PARP1"] <- .095
  }
  gis <- gis[order(gis$Highlight, abs(gis$rho)),]
  
  cr <- cor(gis$InteractionScore.Tau[gis$Highlight], gis$rho[gis$Highlight])
  
  plt <- ggplot(gis, aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
    geom_abline(alpha = .7) +
    geom_abline(slope = 1, intercept = nu.bound, alpha = 0.5, linetype = "dashed") +
    geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 0, alpha = 0.7) +
    geom_hline(yintercept = tau.bound, alpha = 0.5, linetype = "dashed") +
    geom_hline(yintercept = -1*tau.bound, alpha = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0, alpha = 0.7) +
    geom_vline(xintercept = gamma.bound, alpha = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -1*(gamma.bound), alpha= 0.5, linetype= "dashed") +
    geom_point(alpha = .5, size = .85,  col = "#999999", shape = 16) +
    geom_point(data = gis[gis$Highlight & gis$Non != "PARP1",],aes(fill = rho), 
               size = 2, alpha = .7, col = "black", shape = 21) +
    geom_point(data = gis[gis$Highlight & gis$Non == "PARP1",],col = "black", size = 2.2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line()) + 
    coord_fixed() +
    scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
    scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
    scale_fill_gradientn(colors = c("#003C30", "#01665E",  "#35978F", "#80CDC1", 
                                    "#E6E6E6", "#8C510A"),
                         values = c(0, .16, .32, .48, .66, .81, 1),
                         breaks = seq(-.9, .3, .3),
                         labels = c(-.9, -.6, -.3, 0, .3)) +
    xlab("") + ylab("")
  
  plt2 <- ggplot(gis[gis$Highlight,], aes(InteractionScore.Gamma, InteractionScore.Tau)) + 
    geom_abline(alpha = .7) +
    geom_abline(slope = 1, intercept = nu.bound, alpha = 0.65, linetype = "dashed") +
    geom_abline(slope = 1, intercept = -1*nu.bound, alpha = 0.65, linetype = "dashed") +
    geom_hline(yintercept = 0, alpha = 0.7) +
    geom_hline(yintercept = tau.bound, alpha = 0.65, linetype = "dashed") +
    geom_hline(yintercept = -1*tau.bound, alpha = 0.65, linetype = "dashed") +
    geom_vline(xintercept = 0, alpha = 0.7) +
    geom_vline(xintercept = gamma.bound, alpha = 0.65, linetype = "dashed") +
    geom_vline(xintercept = -1*(gamma.bound), alpha= 0.65, linetype= "dashed") +
    #geom_point(alpha = .5, size = .8,  col = "#999999") +
    geom_point(data = gis[gis$Highlight & gis$Non != "PARP1",],aes(fill = rho), 
               size = 2, alpha = .7, col = "black", shape = 21) +
    geom_point(data = gis[gis$Highlight & gis$Non == "PARP1",],col = "black", size = 2.2) +
    annotate(geom = "text", x = 4.5, y = -6, hjust = 0,
             label = paste0("r = ", round(cr, 3))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.length = unit(.2, "cm")) + 
    coord_fixed() +
    scale_x_continuous(limits = c(-14,11), breaks = seq(-12, 9, 3)) +
    scale_y_continuous(limits = c(-19,19), breaks = seq(-18, 18, 3)) +
    scale_fill_gradientn(colors = c("#003C30", "#01665E",  "#35978F", "#80CDC1", 
                                    "#E6E6E6", "#8C510A"),
                         values = c(0, .16, .32, .48, .66, .81, 1),
                         breaks = seq(-.9, .3, .3),
                         labels = c(-.9, -.6, -.3, 0, .3)) +
    xlab("") + ylab("")
  
  
  return(list("plot" = plt, "plotlb" = plt2))
}

source("helper_functions.R")

con_path_pfx <- "../DATA/Interaction_Scores/Compiled/Construct_Scores/"
con.gamma.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_gamma_oi_avg.txt", sep = ""))
con.tau.oi.avg <- read_tsv(paste(con_path_pfx, "all_interaction_scores_tau_oi_avg.txt", sep = ""))
gene.id.map <- read_tsv("../DATA/genecombination_id_map.txt")

rho.single.pheno <- read_tsv("../DATA/single_sgRNA_rho_phenotypes.txt")
rho.single.pheno$gene <- gsub("_.*", "", rho.single.pheno$sgRNA.ID)

rho.single.gene.pheno <- rho.single.pheno %>%
  group_by(gene) %>%
  summarise(rho = mean(Rho.Avg))

# pheno.gene <- pheno %>%
#   group_by(GeneCombinationID) %>%
#   summarise(gamma = mean(Gamma.Avg),
#             tau = mean(Tau.Avg),
#             rho = mean(Rho.Avg))

con.oi.avg <- compute_construct_differential_scores(con.gamma.oi.avg, con.tau.oi.avg)

gamma.gene.gi <- compute_gene_interaction_scores(con.gamma.oi.avg, "GI.z")
tau.gene.gi <- compute_gene_interaction_scores(con.tau.oi.avg, "GI.z")

gamma.bound <- 1.524
gamma.hc.bound <- 2.723
tau.bound <- 1.527
tau.hc.bound <- 2.853
nu.bound <- 2.975

gi.cmp <-inner_join(gamma.gene.gi[,1:3], 
                    tau.gene.gi[,1:3], 
                    by = c("GeneCombinationID", "Category"), suffix = c(".Gamma", ".Tau"))
gi.cmp2 <- inner_join(gi.cmp,
                      gene.id.map) %>%
  filter(Category  == "X+Y")

parp1 <- gtplot("PARP1", gi.cmp2)
dnph1 <- gtplot("DNPH1", gi.cmp2)
brca2 <- gtplot("BRCA2", gi.cmp2)
rif1 <- gtplot("RIF1", gi.cmp2)
fanca <- gtplot("FANCA", gi.cmp2)

parp2 <- gtplot("PARP2", gi.cmp2)
aunip <- gtplot("AUNIP", gi.cmp2)

rnaseh2c <- gtplot("RNASEH2C", gi.cmp2)

babam <- gtplot("BABAM1", gi.cmp2)
uimc <- gtplot("UIMC1", gi.cmp2)
fam <- gtplot("FAM175A", gi.cmp2)
brc <- gtplot("BRCC3", gi.cmp2)

ggsave("../FIGURES/FIGURE-4/parp1_highlight.png", parp1[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/rif1_highlight.png", rif1[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/fanca_highlight.png", fanca[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/dnph1_highlight.png", dnph1[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/brca2_highlight.png", brca2[[1]], width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-S4/parp2_highlight.png", parp2[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-S6/aunip_highlight.png", aunip[[1]], width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-4/parp1_highlight_labels.svg", parp1[[2]], width = 6, height = 10)
ggsave("../FIGURES/FIGURE-4/rif1_highlight_labels.svg", rif1[[2]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/fanca_highlight_labels.svg", fanca[[2]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/dnph1_highlight_labels.svg", dnph1[[2]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-4/brca2_highlight_labels.svg", brca2[[2]], width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-S4/parp2_highlight_labels.svg", parp2[[2]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-S6/aunip_highlight_labels.svg", aunip[[2]], width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/TALK/rnaseh2c_highlight.png", rnaseh2c[[1]], width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-S6/babam1_highlight.png", babam[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-S6/uimc1_highlight.png", uimc[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-S6/fam175a_highlight.png", fam[[1]], width = 6, height = 10, dpi = 300)
ggsave("../FIGURES/FIGURE-S6/brcc3_highlight.png", brc[[1]], width = 6, height = 10, dpi = 300)

ggsave("../FIGURES/FIGURE-S6/babam1_highlight_labels.svg", babam[[2]], width = 6, height = 10)
ggsave("../FIGURES/FIGURE-S6/uimc1_highlight_labels.svg", uimc[[2]], width = 6, height = 10)
ggsave("../FIGURES/FIGURE-S6/fam175a_highlight_labels.svg", fam[[2]], width = 6, height = 10)
ggsave("../FIGURES/FIGURE-S6/brcc3_highlight_labels.svg", brc[[2]], width = 6, height = 10)
