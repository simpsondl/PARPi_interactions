library(ComplexHeatmap)
library(readxl)
library(readr)

horlbeck <- read_excel("F:/Downloads/NIHMS973971-supplement-9.xlsx")
annos <- read_tsv("../DATA/External_Data/screen_annotations_w_gene_inclusion_info.txt")
herken <- read_excel("../../Table S3.xlsx")

colnames(herken)[1] <- "Gene"

horlbeck.genes <- unique(horlbeck$`gene name`)
is.genes <- unique(annos$MapGeneName)
herken.genes <- unique(herken$Gene)

lt <- list(Horlbeck = horlbeck.genes,
           ThisStudy = is.genes)
m1 <- make_comb_mat(lt, mode = "intersect")
UpSet(m1)
