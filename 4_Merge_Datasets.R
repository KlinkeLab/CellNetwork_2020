# This file needs to merge the different digital cytometry measures together
# Create a "Cancer" column with 1 if the sample was cancer and 0 if normal
# create pM0, pM1, and pM2 metrics
# log transform some of the metrics
setwd("~/Documents/Publications/BayesNetwork2/R")
rm(list = ls())

# Read in Proliferation state metric
Prolif_SM <- read.csv(file = "./BRCA_proliferation_Zscore.csv", head=TRUE, stringsAsFactors = FALSE)
Prolif_SM$Cancer <- sapply(Prolif_SM$ID, function(x) ifelse(substr(x,14,15)== "11", 0, 1))

# Read in EMT state metrics
EMT_SM <- read.csv(file = "./BRCA-EMT-SM.csv", head=TRUE, stringsAsFactors = FALSE)
EMT_SM$X <- NULL
DC_TOT <- merge(Prolif_SM, EMT_SM, by.x = "ID", by.y = "TCGAName")

# Read in Digital Cytometry based on LM22 signatures
DC_LM22 <- read.csv(file = "./BRCA_CIBERSORTx_LM22_April2020.csv", head=TRUE, stringsAsFactors = FALSE)
DC_LM22$P.value <- NULL
DC_LM22$Correlation <- NULL
DC_LM22$RMSE<- NULL
DC_LM22$Absolute.score..sig.score. <- NULL

# Find minimum non-zero value
min(DC_LM22$T.cells.CD8[DC_LM22$T.cells.CD8 > 0])
# 0.00336 - so have zero in log transform as 0.003
DC_LM22$T.cells.CD8_lg <- log2(DC_LM22$T.cells.CD8 + 0.003)
plot(density(log2(DC_LM22$T.cells.CD8 + 0.003), adj = 0.5))

# Find minimum non-zero value
min(DC_LM22$NK.cells.activated[DC_LM22$NK.cells.activated > 0])
# 0.001657 - so have zero in log transform as 0.001
DC_LM22$NK.cells.activated_lg <- log2(DC_LM22$NK.cells.activated + 0.001)
plot(density(log2(DC_LM22$NK.cells.activated + 0.001), adj = 0.5))

# Find minimum non-zero value
min(DC_LM22$B.cells.naive[DC_LM22$B.cells.naive > 0])
# 0.000969 - so have zero in log transform as 0.0005
DC_LM22$B.cells.naive_lg <- log2(DC_LM22$B.cells.naive + 0.0009)
plot(density(log2(DC_LM22$B.cells.naive + 0.0009), adj = 0.5))

# Find minimum non-zero value
min(DC_LM22$Neutrophils[DC_LM22$Neutrophils > 0])
# 0.000858 - so have zero in log transform as 0.0005
DC_LM22$Neutrophils_lg <- log2(DC_LM22$Neutrophils + 0.0008)
plot(density(log2(DC_LM22$Neutrophils + 0.0008), adj = 0.5))

# Generate Macrophage orientation measures
DC_LM22$pM0 <- DC_LM22$Macrophages.M0/(DC_LM22$Macrophages.M0 + DC_LM22$Macrophages.M1 + DC_LM22$Macrophages.M2)
DC_LM22$pM1 <- DC_LM22$Macrophages.M1/(DC_LM22$Macrophages.M0 + DC_LM22$Macrophages.M1 + DC_LM22$Macrophages.M2)
DC_LM22$pM2 <- DC_LM22$Macrophages.M2/(DC_LM22$Macrophages.M0 + DC_LM22$Macrophages.M1 + DC_LM22$Macrophages.M2)
DC_TOT <- merge(DC_TOT, DC_LM22, by.x = "ID", by.y = "Mixture")

# Read in Digital Cytometry based on scRNAseq signatures extracted from Tirosh et al. analysis of melanoma tissues
DC_scRNA <- read.csv(file = "./BRCA_CIBERSORTx_scRNA_April2020.csv", head=TRUE, stringsAsFactors = FALSE)
DC_scRNA$P.value <- NULL
DC_scRNA$Correlation <- NULL
DC_scRNA$RMSE<- NULL
DC_scRNA$Absolute.score..sig.score. <- NULL
DC_scRNA$T.cells.CD8.sc <- DC_scRNA$T.cells.CD8
DC_scRNA$T.cells.CD8 <- NULL
# Find minimum non-zero value
min(DC_scRNA$Endothelial.cells[DC_scRNA$Endothelial.cells > 0])
# 0.001146 - so have zero in log transform as 0.001
DC_scRNA$Endothelial.cells_lg <- log2(DC_scRNA$Endothelial.cells + 0.001)
DC_scRNA$CAF_lg <- log2(DC_scRNA$CAF + 0.001)
DC_scRNA$Macrophages_sc_lg <- log2(DC_scRNA$Macrophages + 0.001)
# Find minimum non-zero value
min(DC_scRNA$T.cells.CD4[DC_scRNA$T.cells.CD4 > 0])
# 0.0008447 - so have zero in log transform as 0.0008
DC_scRNA$CD4Tcell_sc_lg <- log2(DC_scRNA$T.cells.CD4 + 0.0008)

plot(density(log2(DC_scRNA$T.cells.CD4 + 0.0008), adj = 0.5))
DC_TOT <- merge(DC_TOT, DC_scRNA, by.x = "ID", by.y = "Mixture")

file <- "BRCA_digitalcytometry_tot"
write.csv(DC_TOT, file = paste("~/Documents/Publications/BayesNetwork2/R/", file, ".csv", sep = ""), row.names = TRUE)



