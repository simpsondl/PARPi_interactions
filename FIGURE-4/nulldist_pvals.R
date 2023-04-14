library(readr)
library(dplyr)
library(ggplot2)
library(qvalue)

# Get file names
fp <- "../DATA/Interaction_Scores/NullDistributionPvalues/gamma"
files <- list.files(fp)
# Last file needs to be imported differently, has 451200 empty lines
last.file <- "allbatches_140002_150001.txt"

for(i in files){
  if(i == last.file){
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE, skip = 558300 )
  } else {
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE)
  }
  
  # These were calculated with wrong total, adjusting
  tmp$X4 <- 2953160 - tmp$X3
  tmp$X6 <- 2953160 - tmp$X5
  
  # Adjust counts where awk returned -1
  tmp$X4[tmp$X3 == -1] <- 0
  tmp$X3[tmp$X3 == -1] <- 2953160
  
  # Fix opposite counts where awk behaved weirdly
  tmp$X6[tmp$X5 == 0 & tmp$X2 < 0] <- 0
  tmp$X5[tmp$X5 == 0 & tmp$X2 < 0] <- 2953160
  
  tmp2 <- tmp %>% 
    group_by(X1) %>% 
    summarise(IS = unique(X2),
              N.Less = sum(X3), 
              N.Greater = sum(X4), 
              N.OppLess = sum(X5), 
              N.OppGreater = sum(X6))
  
  if(i == files[1]){
    all.counts <- tmp2
  } else {
    all.counts <- rbind(all.counts, tmp2)
  }
}



#all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Less + all.counts$N.Greater)
#all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Less + all.counts$N.Greater)
#all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0])/(all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.Greater[all.counts$IS <= 0])
#all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0])/(all.counts$N.Less[all.counts$IS > 0] + all.counts$N.Greater[all.counts$IS > 0])
all.counts$N.Twosided[all.counts$IS <= 0] <- all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0]
all.counts$N.Twosided[all.counts$IS > 0] <- all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0]
all.counts$N.Total <- all.counts$N.Less + all.counts$N.Greater
all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Total)
all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Total)
all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0] + 1)/(all.counts$N.Total[all.counts$IS <= 0] + 1)
all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0] + 1)/(all.counts$N.Total[all.counts$IS > 0] + 1)

qobj <- qvalue(p = all.counts$Ptwosided, fdr.level = .05, lambda = seq(0, .95, .01), pi0.method = "bootstrap")
all.counts$Qvalue <- qobj$qvalues
all.counts$LocalFDR <- qobj$lfdr
all.counts$Significant <- qobj$significant
colnames(all.counts)[1] <- "GeneCombinationID"

all.counts2 <- all.counts[,c("GeneCombinationID", "IS", "N.Twosided", "N.Total", "Ptwosided", "Qvalue", "LocalFDR", "Significant")]
all.counts2 <- all.counts2[order(abs(all.counts2$IS), decreasing = TRUE),]

write_tsv(all.counts2, "../DATA/nulldist_pvals_gamma_genegeneinteractionscores.txt")



# Get file names
fp <- "../DATA/Interaction_Scores/NullDistributionPvalues/tau"
files <- list.files(fp)
# Last file needs to be imported differently, has 342900 empty lines
last.file <- "allbatches_140002_150001.txt"

for(i in files){
  if(i == last.file){
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE, skip = 342900 )
  } else {
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE)
  }
  
  # These were calculated with wrong total, adjusting
  tmp$X4 <- 2953160 - tmp$X3
  tmp$X6 <- 2953160 - tmp$X5
  
  # Adjust counts where awk returned -1
  tmp$X4[tmp$X3 == -1] <- 0
  tmp$X3[tmp$X3 == -1] <- 2953160
  
  # Fix opposite counts where awk behaved weirdly
  tmp$X6[tmp$X5 == 0 & tmp$X2 < 0] <- 0
  tmp$X5[tmp$X5 == 0 & tmp$X2 < 0] <- 2953160
  
  tmp2 <- tmp %>% 
    group_by(X1) %>% 
    summarise(IS = unique(X2),
              N.Less = sum(X3), 
              N.Greater = sum(X4), 
              N.OppLess = sum(X5), 
              N.OppGreater = sum(X6))
  
  if(i == files[1]){
    all.counts <- tmp2
  } else {
    all.counts <- rbind(all.counts, tmp2)
  }
}

all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0])/(all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.Greater[all.counts$IS <= 0])
all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0])/(all.counts$N.Less[all.counts$IS > 0] + all.counts$N.Greater[all.counts$IS > 0])
all.counts$N.Twosided[all.counts$IS <= 0] <- all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0]
all.counts$N.Twosided[all.counts$IS > 0] <- all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0]
all.counts$N.Total <- all.counts$N.Less + all.counts$N.Greater

qobj <- qvalue(p = all.counts$Ptwosided, fdr.level = .05, lambda = seq(0, .95, .01), pi0.method = "bootstrap")
all.counts$Qvalue <- qobj$qvalues
all.counts$LocalFDR <- qobj$lfdr
all.counts$Significant <- qobj$significant

all.counts2 <- all.counts[,c(1:4,10:11,7:9,12:13)]
all.counts2 <- all.counts2[order(abs(all.counts2$IS), decreasing = TRUE),]

write_tsv(all.counts2, "../DATA/nulldist_pvals_tau_genegeneinteractionscores.txt")




# Get file names
fp <- "../DATA/Interaction_Scores/NullDistributionPvalues/nu"
files <- list.files(fp)
# Last file needs to be imported differently, has 612100 empty lines
last.file <- "allbatches_140002_150001.txt"

for(i in files){
  if(i == last.file){
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE, skip = 612100 )
  } else {
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE)
  }
  
  # These were calculated with wrong total, adjusting
  tmp$X4 <- 2877580 - tmp$X3
  tmp$X6 <- 2877580 - tmp$X5
  
  # Adjust counts where awk returned -1
  tmp$X4[tmp$X3 == -1] <- 0
  tmp$X3[tmp$X3 == -1] <- 2877580
  
  # Fix opposite counts where awk behaved weirdly
  tmp$X6[tmp$X5 == 0 & tmp$X2 < 0] <- 0
  tmp$X5[tmp$X5 == 0 & tmp$X2 < 0] <- 2877580
  
  tmp2 <- tmp %>% 
    group_by(X1) %>% 
    summarise(IS = unique(X2),
              N.Less = sum(X3), 
              N.Greater = sum(X4), 
              N.OppLess = sum(X5), 
              N.OppGreater = sum(X6))
  
  if(i == files[1]){
    all.counts <- tmp2
  } else {
    all.counts <- rbind(all.counts, tmp2)
  }
}

all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0])/(all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.Greater[all.counts$IS <= 0])
all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0])/(all.counts$N.Less[all.counts$IS > 0] + all.counts$N.Greater[all.counts$IS > 0])
all.counts$N.Twosided[all.counts$IS <= 0] <- all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0]
all.counts$N.Twosided[all.counts$IS > 0] <- all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0]
all.counts$N.Total <- all.counts$N.Less + all.counts$N.Greater

qobj <- qvalue(p = all.counts$Ptwosided, fdr.level = .05, lambda = seq(0, .95, .01), pi0.method = "bootstrap")
all.counts$Qvalue <- qobj$qvalues
all.counts$LocalFDR <- qobj$lfdr
all.counts$Significant <- qobj$significant

all.counts2 <- all.counts[,c(1:4,10:11,7:9,12:13)]
all.counts2 <- all.counts2[order(abs(all.counts2$IS), decreasing = TRUE),]

write_tsv(all.counts2, "../DATA/nulldist_pvals_nu_genegeneinteractionscores.txt")


######################################################################

# Get file names
fp <- "../DATA/Interaction_Scores/NullDistributionPvalues/gamma_r1/"
files <- list.files(fp)
# Last file needs to be imported differently, has 342900 empty lines
last.file <- "allbatches_140002_150001.txt"

for(i in files){
  if(i == last.file){
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE, skip = 1116600 )
  } else {
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE)
  }
  
  # These were calculated with wrong total, adjusting
  #tmp$X4 <- 2953160 - tmp$X3
  #tmp$X6 <- 2953160 - tmp$X5
  
  # Adjust counts where awk returned -1
  tmp$X4[tmp$X3 == -1] <- 0
  tmp$X3[tmp$X3 == -1] <- 147658
  
  # Fix opposite counts where awk behaved weirdly
  tmp$X6[tmp$X5 == 0 & tmp$X2 < 0] <- 0
  tmp$X5[tmp$X5 == 0 & tmp$X2 < 0] <- 147658
  
  tmp2 <- tmp %>% 
    group_by(X1) %>% 
    summarise(IS = unique(X2),
              N.Less = sum(X3), 
              N.Greater = sum(X4), 
              N.OppLess = sum(X5), 
              N.OppGreater = sum(X6))
  
  if(i == files[1]){
    all.counts <- tmp2
  } else {
    all.counts <- rbind(all.counts, tmp2)
  }
}

all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0])/(all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.Greater[all.counts$IS <= 0])
all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0])/(all.counts$N.Less[all.counts$IS > 0] + all.counts$N.Greater[all.counts$IS > 0])
all.counts$N.Twosided[all.counts$IS <= 0] <- all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0]
all.counts$N.Twosided[all.counts$IS > 0] <- all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0]
all.counts$N.Total <- all.counts$N.Less + all.counts$N.Greater

qobj <- qvalue(p = all.counts$Ptwosided, fdr.level = .05, lambda = seq(0, .95, .01), pi0.method = "bootstrap")
all.counts$Qvalue <- qobj$qvalues
all.counts$LocalFDR <- qobj$lfdr
all.counts$Significant <- qobj$significant

all.counts2 <- all.counts[,c(1:4,10:11,7:9,12:13)]
all.counts2 <- all.counts2[order(abs(all.counts2$IS), decreasing = TRUE),]

write_tsv(all.counts2, "../DATA/Interaction_Scores/NullDistributionPvalues/gamma_r1/gamma_r1_qvalues.txt")

######################################################################

# Get file names
fp <- "../DATA/Interaction_Scores/NullDistributionPvalues/gamma_r2/"
files <- list.files(fp)
# Last file needs to be imported differently, has 342900 empty lines
last.file <- "allbatches_140002_150001.txt"

for(i in files){
  if(i == last.file){
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE, skip = 1116600 )
  } else {
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE)
  }
  
  # These were calculated with wrong total, adjusting
  #tmp$X4 <- 2953160 - tmp$X3
  #tmp$X6 <- 2953160 - tmp$X5
  
  # Adjust counts where awk returned -1
  tmp$X4[tmp$X3 == -1] <- 0
  tmp$X3[tmp$X3 == -1] <- 147658
  
  # Fix opposite counts where awk behaved weirdly
  tmp$X6[tmp$X5 == 0 & tmp$X2 < 0] <- 0
  tmp$X5[tmp$X5 == 0 & tmp$X2 < 0] <- 147658
  
  tmp2 <- tmp %>% 
    group_by(X1) %>% 
    summarise(IS = unique(X2),
              N.Less = sum(X3), 
              N.Greater = sum(X4), 
              N.OppLess = sum(X5), 
              N.OppGreater = sum(X6))
  
  if(i == files[1]){
    all.counts <- tmp2
  } else {
    all.counts <- rbind(all.counts, tmp2)
  }
}

all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0])/(all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.Greater[all.counts$IS <= 0])
all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0])/(all.counts$N.Less[all.counts$IS > 0] + all.counts$N.Greater[all.counts$IS > 0])
all.counts$N.Twosided[all.counts$IS <= 0] <- all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0]
all.counts$N.Twosided[all.counts$IS > 0] <- all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0]
all.counts$N.Total <- all.counts$N.Less + all.counts$N.Greater

qobj <- qvalue(p = all.counts$Ptwosided, fdr.level = .05, lambda = seq(0, .95, .01), pi0.method = "bootstrap")
all.counts$Qvalue <- qobj$qvalues
all.counts$LocalFDR <- qobj$lfdr
all.counts$Significant <- qobj$significant

all.counts2 <- all.counts[,c(1:4,10:11,7:9,12:13)]
all.counts2 <- all.counts2[order(abs(all.counts2$IS), decreasing = TRUE),]

write_tsv(all.counts2, "../DATA/Interaction_Scores/NullDistributionPvalues/gamma_r2/gamma_r2_qvalues.txt")


######################################################################

# Get file names
fp <- "../DATA/Interaction_Scores/NullDistributionPvalues/gamma_avg/"
files <- list.files(fp)
# Last file needs to be imported differently, has 342900 empty lines
last.file <- "allbatches_140002_150001.txt"

for(i in files){
  if(i == last.file){
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE, skip = 1116600 )
  } else {
    tmp <- read_table(paste(fp, i, sep = "/"), col_names = FALSE)
  }
  
  # These were calculated with wrong total, adjusting
  #tmp$X4 <- 2953160 - tmp$X3
  #tmp$X6 <- 2953160 - tmp$X5
  
  # Adjust counts where awk returned -1
  tmp$X4[tmp$X3 == -1] <- 0
  tmp$X3[tmp$X3 == -1] <- 147658
  
  # Fix opposite counts where awk behaved weirdly
  tmp$X6[tmp$X5 == 0 & tmp$X2 < 0] <- 0
  tmp$X5[tmp$X5 == 0 & tmp$X2 < 0] <- 147658
  
  tmp2 <- tmp %>% 
    group_by(X1) %>% 
    summarise(IS = unique(X2),
              N.Less = sum(X3), 
              N.Greater = sum(X4), 
              N.OppLess = sum(X5), 
              N.OppGreater = sum(X6))
  
  if(i == files[1]){
    all.counts <- tmp2
  } else {
    all.counts <- rbind(all.counts, tmp2)
  }
}

all.counts$Pless <- (all.counts$N.Less)/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Pgreater <- all.counts$N.Greater/(all.counts$N.Less + all.counts$N.Greater)
all.counts$Ptwosided[all.counts$IS <= 0] <- (all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0])/(all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.Greater[all.counts$IS <= 0])
all.counts$Ptwosided[all.counts$IS > 0] <- (all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0])/(all.counts$N.Less[all.counts$IS > 0] + all.counts$N.Greater[all.counts$IS > 0])
all.counts$N.Twosided[all.counts$IS <= 0] <- all.counts$N.Less[all.counts$IS <= 0] + all.counts$N.OppGreater[all.counts$IS <= 0]
all.counts$N.Twosided[all.counts$IS > 0] <- all.counts$N.Greater[all.counts$IS > 0] + all.counts$N.OppLess[all.counts$IS > 0]
all.counts$N.Total <- all.counts$N.Less + all.counts$N.Greater

qobj <- qvalue(p = all.counts$Ptwosided, fdr.level = .05, lambda = seq(0, .95, .01), pi0.method = "bootstrap")
all.counts$Qvalue <- qobj$qvalues
all.counts$LocalFDR <- qobj$lfdr
all.counts$Significant <- qobj$significant

all.counts2 <- all.counts[,c(1:4,10:11,7:9,12:13)]
all.counts2 <- all.counts2[order(abs(all.counts2$IS), decreasing = TRUE),]

write_tsv(all.counts2, "../DATA/Interaction_Scores/NullDistributionPvalues/gamma_avg/gamma_avg_qvalues.txt")
