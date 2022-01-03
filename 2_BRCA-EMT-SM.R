# Calculate values for Epithelial and Mesenchymal state metrics
setwd("E:/Lab Files/Bayes Network/CellNetwork_2020-master")
rm(list = ls())

load(file = "./BRCA_TPM_HK.Rda")

EMTsig <- read.table(file = "./EMT_BRCA.csv", sep = ",", head=TRUE, colClasses = c("character", "numeric", "factor"))

#Epithelial signature
# Changing the following gene symbols in the .csv file may be needed -> SEPP1:SELENOP, C1orf106:INAVA

Esig <- EMTsig[EMTsig$State == "E",]

colnames(Esig) <- c('GENE_SYMBOL', 'Ki_log2_TPM', 'State' )

RNAseq.E <- TCGA.RNAseq.HK[match(Esig$GENE_SYMBOL, row.names(TCGA.RNAseq.HK)), ]

E.thresh <- kmeans(t(log2(RNAseq.E + 0.03)), centers = 2)$centers
E.thresh2 <- apply(E.thresh, 2, mean)

#Check to see that Ki's are the same
plot(E.thresh2, Esig$Ki_log2_TPM)
lines(c(-2,8), c(-2,8))

EBrCa <- rep(0, dim(RNAseq.E)[2])
for (i in 1:dim(RNAseq.E)[2])
{
  EBrCa[i] <- mean(RNAseq.E[,i]/(RNAseq.E[,i] + 2^E.thresh2), na.rm = TRUE)
}

# Mesenchymal signature
# Changing the following gene symbols in the .csv file may be needed -> C7orf10:SUGCT, LEPRE1:P3H1, LHFP:LHFPL6, WISP1:CCN4

Msig <- EMTsig[EMTsig$State == "M",]
colnames(Msig) <- c('GENE_SYMBOL', 'Ki_log2_TPM', 'State' )

RNAseq.M <- TCGA.RNAseq.HK[match(Msig$GENE_SYMBOL, row.names(TCGA.RNAseq.HK)), ]

M.thresh <- kmeans(t(log2(RNAseq.M + 0.03)), centers = 2)$centers
M.thresh2 <- apply(M.thresh, 2, mean)

#Check to see that Ki's are the same
plot(M.thresh2, Msig$Ki_log2_TPM)
lines(c(-2,8), c(-2,8))

MBrCa <- rep(0, dim(RNAseq.M)[2])
for (i in 1:dim(RNAseq.M)[2])
{
  MBrCa[i] <- mean(RNAseq.M[,i]/(RNAseq.M[,i] + 2^M.thresh2), na.rm = TRUE)
}

#Combine Epithelial and Mesenchymal genes for each sample
# If you changed WISP1 to CCN4 in the .csv, "WISP1" to "CCN4" in line 51 as well.
BrCa_State <- data.frame(TCGAName = colnames(RNAseq.M), Epithelial = EBrCa, Mesenchymal = MBrCa, CCN4 = log2(as.numeric(RNAseq.M["CCN4",]) + 0.001))
write.csv(BrCa_State, file = "./BRCA-EMT-SM.csv")
