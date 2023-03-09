library(dplyr)
library(tidyverse)
library(strataG)
source("R/gelato.KKM.R")

load("data/hap.freqs.rda")
hap.freqs$CentAm <- hap.freqs$CentAm + hap.freqs$GRO.OAX

MCP.by.year.hapfreq <- read.csv("data-raw/MCP-by-year-hapfreqs.csv", header = TRUE, check.names = FALSE)
yrs <- names(MCP.by.year.hapfreq)[-1]
hap.freqs <- left_join(hap.freqs[,1:4], MCP.by.year.hapfreq)
hap.freqs[is.na(hap.freqs)] <- 0

# Convert hap.freqs into dataframe strataG can use (row=ind, cols=id, hap, stratum)
hap.df <- do.call('rbind', lapply(2:ncol(hap.freqs), function(i){
  inds <- data.frame(do.call('c', lapply(1:nrow(hap.freqs), function(h){
    rep(hap.freqs$Haplotype[h], hap.freqs[h,i])
  })))
  cbind(names(hap.freqs)[i], inds)
}))

hap.df <- cbind(paste0("ind.", seq(1:nrow(hap.df))), hap.df)
names(hap.df) <- c("id","stratum","hap")

# Read in haplotype sequences
seqs <- read.fasta("data-raw/MCP-new-aligned-to-Baker.fasta")

# Create gtypes object
hap.g <- df2gtypes(hap.df, ploidy = 1, strata.col = 2, loc.col = 3, sequences = seqs)

# GELATo analysis
gelato.by.yr <- lapply(yrs, function(y){
  g.y <- df2gtypes(filter(hap.df, stratum %in% c("MMex", "CentAm", y)), ploidy = 1, strata.col = 2, loc.col = 3, sequences = seqs)
  hap.gel <- gelato.KKM(g.y, unknown.strata = y, nrep = 1000)
})

# Format and save results
gelato.ass.prob <- do.call(rbind, lapply(gelato.by.yr, function(y){
  y$assign.prob
}))
gelato.ass.prob$odds <- gelato.ass.prob$MMex/gelato.ass.prob$CentAm

gelato.lnL <- do.call(rbind, lapply(gelato.by.yr, function(y){
  data.frame(cbind(MMex = y$likelihoods[[1]]$MMex$log.Lik.smry[3],
             CentAm = y$likelihoods[[1]]$CentAm$log.Lik.smry[3]))
}))
rownames(gelato.lnL) <- yrs

save(gelato.ass.prob, gelato.by.yr, gelato.lnL, file = "results/gelato.and.fst.by.year.rda")

write.csv(gelato.lnL, file = "results/gelato.lnL.csv")
write.csv(gelato.ass.prob, file = "results/gelato.ass.prob.csv")

