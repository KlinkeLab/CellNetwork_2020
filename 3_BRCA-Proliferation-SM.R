# 3_BRCA-Proliferation-SM.R - calculate proliferation state metric
# Copyright (C) 2022 Klinke Lab
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(colorspace)
library(MASS)
setwd("E:/Lab Files/Bayes Network/CellNetwork_2020-master")
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
