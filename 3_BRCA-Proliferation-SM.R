library(colorspace)
library(MASS)
setwd("~/Documents/Publications/BayesNetwork2/R")
rm(list = ls())

# Get BRCA sample IDs
load(file = "./BRCA_TPM_HK.Rda")

#pos.list
Glist <- read.csv("./proliferation_BRCA.csv", head = TRUE, na.strings = "NA", colClasses = c("character"))

ar.pro <- TCGA.RNAseq.HK[match(Glist$Gene_Symbol, rownames(TCGA.RNAseq.HK)),]

summary(Glist$Gene_Symbol %in% rownames(TCGA.RNAseq.HK))

ar.pos <- t(log2(as.matrix(ar.pro)+0.001))
#ar.pos <- t(as.matrix(ar.pro))

pMPC <- apply(ar.pos, 2, median)
pSPC <- apply(ar.pos, 2, sd)

# modified way
ar.d <- data.frame(ID = rownames(ar.pos), proliferation = 0, stringsAsFactors = FALSE)
factor = 3

for(a in 1:ncol(ar.pos))
{
  for(b in 1:nrow(ar.d))
  {
    ar.d$proliferation[b] = ar.d$proliferation[b] + (ar.pos[b,a] - pMPC[a])/pSPC[a] # Z-score
  }
}

write.csv(ar.d, paste("./BRCA_proliferation_Zscore.csv", sep = ""), row.names=FALSE)

# Plot distributions in proliferation separately for tumor and normal groups
N_DC <- ar.d$proliferation[substr(ar.d$ID,14,15)== "11"]
T_DC <- ar.d$proliferation[substr(ar.d$ID,14,15)== "01"]

Nden = density(N_DC, adj = 0.2)
Tden = density(T_DC, adj = 0.2)
xlimit = c(min(ar.d$proliferation), max(ar.d$proliferation))

plot(Nden$x, Nden$y, type = "l", xlim = xlimit, xlab = "Metric", ylab = "Density", col = "black", main = list(cex = 0.8, "Proliferation"))
lines(Tden$x, Tden$y, col = "red")
