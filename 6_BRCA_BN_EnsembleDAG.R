library(colorspace)
library(stats) 
library(bnlearn)
#library(lattice) 
#library(Rgraphviz)

#library(MASS)
setwd("~/Documents/Publications/BayesNetwork2/R")
rm(list = ls())

# Making compute cluster
cl = makeCluster(6, type = "SOCK")
rand = clusterEvalQ(cl, runif(10))

# Taking ensemble DAG
file <- "BRCA_digitalcytometry_tot_noise"
test <- read.csv(paste("./", file, ".csv", sep = ""), head = TRUE)

DM <- c("Cancer", "Endothelial.cells_lg", "CAF_lg", "T.cells.CD8_lg", "NK.cells_lg", "Macrophages_sc_lg", 
        "CD4Tcell_sc_lg", "B.cells.naive_lg", "Neutrophils_lg", "proliferation", "Epithelial", "Mesenchymal", 
        "CCN4", "pM0", "pM1", "pM2")

test2 <- test[, DM]
test2$Cancer <- as.numeric(test$Cancer)

test2 <- test2[complete.cases(test2),]

tmp <- discretize(test2, method = "interval", breaks = c(2,rep(15,15)))
dtest2 <- as.data.frame(lapply(tmp, function(x) (as.numeric(x) - 1)/max(as.numeric(x) -1)))

BRCA2 <- dtest2

corrcoef <- cor(dtest2)
colnames(corrcoef) <- colnames(dtest2)
rownames(corrcoef) <- colnames(dtest2)

# Only arcs into CD8 T cells (leaf node), only arcs out from Cancer (root node), mostly arcs out of CCN4 (with exception for Cancer), 
# Only arcs into CD4 T cells as number of zeros is high
# Only arcs into Neutrophils as number of zeros is high

CCN4black = data.frame(from = c("Endothelial.cells_lg", "CAF_lg", "T.cells.CD8_lg", "NK.cells_lg", "Macrophages_sc_lg", 
                                "CD4Tcell_sc_lg", "B.cells.naive_lg", "Neutrophils_lg", "proliferation", "Epithelial", "Mesenchymal", 
                                "CCN4", "pM0", "pM1", "pM2", 
                                "Endothelial.cells_lg", "CAF_lg", "T.cells.CD8_lg", "NK.cells_lg", "Macrophages_sc_lg", 
                                "CD4Tcell_sc_lg", "B.cells.naive_lg", "Neutrophils_lg", "proliferation", "Epithelial", "Mesenchymal", 
                                "pM0", "pM1", "pM2", "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg",
                                "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg",
                                "T.cells.CD8_lg", "T.cells.CD8_lg", "T.cells.CD8_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", 
                                "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", 
                                "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", 
                                "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg", 
                                "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg", 
                                "Neutrophils_lg", "Neutrophils_lg", "Neutrophils_lg"), 
                       to = c("Cancer", "Cancer", "Cancer", "Cancer", "Cancer", 
                              "Cancer", "Cancer", "Cancer", "Cancer", "Cancer", 
                              "Cancer", "Cancer", "Cancer", "Cancer", "Cancer",  
                              "CCN4", "CCN4", "CCN4", "CCN4", "CCN4", 
                              "CCN4", "CCN4", "CCN4", "CCN4", "CCN4", 
                              "CCN4", "CCN4", "CCN4", "CCN4", "Endothelial.cells_lg", "CAF_lg", "NK.cells_lg", "Macrophages_sc_lg", 
                              "CD4Tcell_sc_lg", "B.cells.naive_lg", "Neutrophils_lg", "proliferation", "Epithelial", "Mesenchymal", 
                              "pM0", "pM1", "pM2", "Endothelial.cells_lg", "CAF_lg", "NK.cells_lg", "Macrophages_sc_lg", 
                              "T.cells.CD8_lg", "B.cells.naive_lg", "Neutrophils_lg", "proliferation", "Epithelial", "Mesenchymal", 
                              "pM0", "pM1", "pM2", "Endothelial.cells_lg", "CAF_lg", "NK.cells_lg", "Macrophages_sc_lg", 
                              "T.cells.CD8_lg", "B.cells.naive_lg", "CD4Tcell_sc_lg", "proliferation", "Epithelial", "Mesenchymal", 
                              "pM0", "pM1", "pM2"))

CCN4white = data.frame(from = c("CCN4", "pM0", "pM0", "Mesenchymal", "Cancer", 
                                "B.cells.naive_lg", "Cancer", "pM1", "Macrophages_sc_lg", "Macrophages_sc_lg", 
                                "CAF_lg", "Cancer", "Cancer", "Endothelial.cells_lg", "B.cells.naive_lg", 
                                "Cancer", "pM1", "Cancer", "pM2", "Cancer", 
                                "CCN4", "Cancer", "CCN4"),
                          to = c("Mesenchymal", "pM2", "pM1", "CAF_lg", "CCN4", 
                                 "T.cells.CD8_lg", "Epithelial", "T.cells.CD8_lg", "T.cells.CD8_lg", "CD4Tcell_sc_lg", 
                                 "T.cells.CD8_lg", "Endothelial.cells_lg", "proliferation", "CD4Tcell_sc_lg", "CD4Tcell_sc_lg", 
                                 "pM2", "CAF_lg", "Mesenchymal", "Neutrophils_lg", "pM1", 
                                 "NK.cells_lg", "B.cells.naive_lg", "Macrophages_sc_lg"))


nboot = 10000 
arstr <- boot.strength(dtest2, R = nboot, algorithm = "mmhc", cluster = cl, algorithm.args = list(whitelist = CCN4white, blacklist = CCN4black))

an.BRCA <- averaged.network(arstr)
ars.BRCA <- arcs(an.BRCA)
CorSign <- rep("+", nrow(ars.BRCA))

fBRCABN <- bn.fit(an.BRCA, data = dtest2)

for(b in 1:nrow(ars.BRCA))
{
  CorSign[b] <- ifelse(corrcoef[match(ars.BRCA[b,1], colnames(corrcoef)), match(ars.BRCA[b,2], colnames(corrcoef))] >0, "+", "-")
}

HLarcs <- ars.BRCA[CorSign == "-",]
BRCA_str = arc.strength(an.BRCA, data = dtest2)

strength.plot(an.BRCA, BRCA_str, shape = "ellipse", highlight = list(arcs = HLarcs))

#Need to include whether the particular edge included in the network has a negative or a positive correlation
# While this is correct normally, ultimately you should look at the sign of the appropriate term in the linear model (fBRCABN)
write.csv(cbind(BRCA_str, CorSign), paste(file, "Bootstrap", nboot, "_mmhc_Hybrid_Merge.csv", sep = ' '), row.names=FALSE)

# stop the cluster when we are done
stopCluster(cl)

pdf("BRCA-ModelvsData-hybrid-April20.pdf", width = 11, height = 6)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
#Comparison plots between BN model and data
# BRCA CCN4 vs CD8 T cells
sim = cpdist(fBRCABN, nodes = c("CCN4", "T.cells.CD8_lg"), n = 10^5, evidence = (Cancer > 0.85))
ColSim <- densCols(x = sim$CCN4, y = sim$T.cells.CD8, nbin = 256, colramp = colorRampPalette(blues9[-(1:2)]))

plot(sim$CCN4, sim$T.cells.CD8, col = ColSim, xlab = "Normalized CCN4 (log2)", ylab = "CD8 T cells metric", xlim = c(0, 1.2), ylim = c(0.0, 1.0))
simData = data.frame(x = sim$CCN4, y = sim$T.cells.CD8)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

coef(lm(y~x, data=simData))

sim2 = cpdist(fBRCABN, nodes = c("CCN4", "T.cells.CD8_lg"), n = 10^6, evidence = (Cancer < 0.15))
Lab.palette <- colorRampPalette(c("red", "orange", "yellow"), space = "Lab")
ColSim2 <- densCols(x = sim2$CCN4, y = sim2$T.cells.CD8, nbin = 256, colramp = Lab.palette)
points(sim2$CCN4, sim2$T.cells.CD8, col = ColSim2)
simData2 = data.frame(x = sim2$CCN4, y = sim2$T.cells.CD8)
abline(coef(lm(y ~ x, data = simData2)), col = "red", lwd = 3)
abline(coef(lm(y ~ x, data = simData2)), col = "orange", lwd = 2)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

minX <- min(test2$CCN4)
maxX <- max(test2$CCN4)
minY <- min(test2$T.cells.CD8)
maxY <- max(test2$T.cells.CD8)

points( (test2$CCN4[test2$Cancer == 1] - minX)/(maxX - minX), (test2$T.cells.CD8[test2$Cancer == 1] - minY)/(maxY - minY), col = "black")
points( (test2$CCN4[test2$Cancer == 0] - minX)/(maxX - minX), (test2$T.cells.CD8[test2$Cancer == 0] - minY)/(maxY - minY), pch = 19, col = "black")


########################

# BRCA CCN4 vs CD4 T cells
sim = cpdist(fBRCABN, nodes = c("CCN4", "CD4Tcell_sc_lg"), n = 10^4, evidence = (Cancer > 0.85))
ColSim <- densCols(x = sim$CCN4, y = sim$CD4Tcell_sc_lg, nbin = 256, colramp = colorRampPalette(blues9[-(1:2)]))

plot(sim$CCN4, sim$CD4Tcell_sc_lg, col = ColSim, xlab = "Normalized CCN4", ylab = "CD4 T cells metric", xlim = c(0, 1.2), ylim = c(0.0, 1.0))
simData = data.frame(x = sim$CCN4, y = sim$CD4Tcell_sc_lg)

coef(lm(y~x, data=simData))

sim2 = cpdist(fBRCABN, nodes = c("CCN4", "CD4Tcell_sc_lg"), n = 10^6, evidence = (Cancer < 0.15))
Lab.palette <- colorRampPalette(c("red", "orange", "yellow"), space = "Lab")
ColSim2 <- densCols(x = sim2$CCN4, y = sim2$CD4Tcell_sc_lg, nbin = 256, colramp = Lab.palette)
points(sim2$CCN4, sim2$CD4Tcell_sc_lg, col = ColSim2)
simData2 = data.frame(x = sim2$CCN4, y = sim2$CD4Tcell_sc_lg)
abline(coef(lm(y ~ x, data = simData2)), col = "red", lwd = 3)
abline(coef(lm(y ~ x, data = simData2)), col = "orange", lwd = 2)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

minX <- min(test2$CCN4)
maxX <- max(test2$CCN4)
minY <- min(test2$CD4Tcell_sc_lg)
maxY <- max(test2$CD4Tcell_sc_lg)

points( (test2$CCN4[test2$Cancer == 1] - minX)/(maxX - minX), (test2$CD4Tcell_sc_lg[test2$Cancer == 1] - minY)/(maxY - minY), col = "black")
points( (test2$CCN4[test2$Cancer == 0] - minX)/(maxX - minX), (test2$CD4Tcell_sc_lg[test2$Cancer == 0] - minY)/(maxY - minY), pch = 19, col = "black")

#################

# BRCA CCN4 vs NK cells
sim = cpdist(fBRCABN, nodes = c("CCN4", "NK.cells_lg"), n = 10^5, evidence = (Cancer > 0.85))
ColSim <- densCols(x = sim$CCN4, y = sim$NK.cells_lg, nbin = 256, colramp = colorRampPalette(blues9[-(1:2)]))

plot(sim$CCN4, sim$NK.cells_lg, col = ColSim, xlab = "Normalized CCN4", ylab = "Active NK cells metric", xlim = c(0, 1.2), ylim = c(0.0, 1.0))
simData = data.frame(x = sim$CCN4, y = sim$NK.cells_lg)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

coef(lm(y~x, data=simData))

sim2 = cpdist(fBRCABN, nodes = c("CCN4", "NK.cells_lg"), n = 10^6, evidence = (Cancer < 0.15))
Lab.palette <- colorRampPalette(c("red", "orange", "yellow"), space = "Lab")
ColSim2 <- densCols(x = sim2$CCN4, y = sim2$NK.cells_lg, nbin = 256, colramp = Lab.palette)
points(sim2$CCN4, sim2$NK.cells_lg, col = ColSim2)
simData2 = data.frame(x = sim2$CCN4, y = sim2$NK.cells_lg)
abline(coef(lm(y ~ x, data = simData2)), col = "red", lwd = 3)
abline(coef(lm(y ~ x, data = simData2)), col = "orange", lwd = 2)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

minX <- min(test2$CCN4)
maxX <- max(test2$CCN4)
minY <- min(test2$NK.cells_lg)
maxY <- max(test2$NK.cells_lg)

points( (test2$CCN4[test2$Cancer == 1] - minX)/(maxX - minX), (test2$NK.cells_lg[test2$Cancer == 1] - minY)/(maxY - minY), col = "black")
points( (test2$CCN4[test2$Cancer == 0] - minX)/(maxX - minX), (test2$NK.cells_lg[test2$Cancer == 0] - minY)/(maxY - minY), pch = 19, col = "black")

###############################################################

# GEO CCN4 vs CAFs
sim = cpdist(fBRCABN, nodes = c("CCN4", "CAF_lg"), n = 10^5, evidence = (Cancer > 0.85))
ColSim <- densCols(x = sim$CCN4, y = sim$CAF_lg, nbin = 256, colramp = colorRampPalette(blues9[-(1:2)]))

plot(sim$CCN4, sim$CAF_lg, col = ColSim, xlab = "Normalized CCN4", ylab = "Cancer Associated Fibroblasts", xlim = c(0, 1.2), ylim = c(0.0, 1.0))
simData = data.frame(x = sim$CCN4, y = sim$CAF_lg)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

coef(lm(y~x, data=simData))

sim2 = cpdist(fBRCABN, nodes = c("CCN4", "CAF_lg"), n = 10^6, evidence = (Cancer < 0.15))
Lab.palette <- colorRampPalette(c("red", "orange", "yellow"), space = "Lab")
ColSim2 <- densCols(x = sim2$CCN4, y = sim2$CAF_lg, nbin = 256, colramp = Lab.palette)
points(sim2$CCN4, sim2$CAF_lg, col = ColSim2)
simData2 = data.frame(x = sim2$CCN4, y = sim2$CAF_lg)
abline(coef(lm(y ~ x, data = simData2)), col = "red", lwd = 3)
abline(coef(lm(y ~ x, data = simData2)), col = "orange", lwd = 2)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

minX <- min(test2$CCN4)
maxX <- max(test2$CCN4)
minY <- min(test2$CAF_lg)
maxY <- max(test2$CAF_lg)

points( (test2$CCN4[test2$Cancer == 1] - minX)/(maxX - minX), (test2$CAF_lg[test2$Cancer == 1] - minY)/(maxY - minY), col = "black")
points( (test2$CCN4[test2$Cancer == 0] - minX)/(maxX - minX), (test2$CAF_lg[test2$Cancer == 0] - minY)/(maxY - minY), pch = 19, col = "black")


###############################################################

# BRCA CCN4 vs B cells
sim = cpdist(fBRCABN, nodes = c("CCN4", "B.cells.naive_lg"), n = 10^5, evidence = (Cancer > 0.85))
ColSim <- densCols(x = sim$CCN4, y = sim$B.cells.naive_lg, nbin = 256, colramp = colorRampPalette(blues9[-(1:2)]))

plot(sim$CCN4, sim$B.cells.naive_lg, col = ColSim, xlab = "Normalized CCN4", ylab = "B cells", xlim = c(0, 1.2), ylim = c(0.0, 1.0))
simData = data.frame(x = sim$CCN4, y = sim$B.cells.naive_lg)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

coef(lm(y~x, data=simData))

sim2 = cpdist(fBRCABN, nodes = c("CCN4", "B.cells.naive_lg"), n = 10^6, evidence = (Cancer < 0.15))
Lab.palette <- colorRampPalette(c("red", "orange", "yellow"), space = "Lab")
ColSim2 <- densCols(x = sim2$CCN4, y = sim2$B.cells.naive_lg, nbin = 256, colramp = Lab.palette)
points(sim2$CCN4, sim2$B.cells.naive_lg, col = ColSim2)
simData2 = data.frame(x = sim2$CCN4, y = sim2$B.cells.naive_lg)
abline(coef(lm(y ~ x, data = simData2)), col = "red", lwd = 3)
abline(coef(lm(y ~ x, data = simData2)), col = "orange", lwd = 2)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

minX <- min(test2$CCN4)
maxX <- max(test2$CCN4)
minY <- min(test2$B.cells.naive_lg)
maxY <- max(test2$B.cells.naive_lg)

points( (test2$CCN4[test2$Cancer == 1] - minX)/(maxX - minX), (test2$B.cells.naive_lg[test2$Cancer == 1] - minY)/(maxY - minY), col = "black")
points( (test2$CCN4[test2$Cancer == 0] - minX)/(maxX - minX), (test2$B.cells.naive_lg[test2$Cancer == 0] - minY)/(maxY - minY), pch = 19, col = "black")


###############################################################

# BRCA CCN4 vs Macrophages
sim = cpdist(fBRCABN, nodes = c("CCN4", "Macrophages_sc_lg"), n = 10^5, evidence = (Cancer > 0.85))
ColSim <- densCols(x = sim$CCN4, y = sim$Macrophages_sc_lg, nbin = 256, colramp = colorRampPalette(blues9[-(1:2)]))

plot(sim$CCN4, sim$Macrophages_sc_lg, col = ColSim, xlab = "Normalized CCN4", ylab = "Macrophages", xlim = c(0, 1.2), ylim = c(0.0, 1.0))
simData = data.frame(x = sim$CCN4, y = sim$Macrophages_sc_lg)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

coef(lm(y~x, data=simData))

sim2 = cpdist(fBRCABN, nodes = c("CCN4", "Macrophages_sc_lg"), n = 10^6, evidence = (Cancer < 0.15))
Lab.palette <- colorRampPalette(c("red", "orange", "yellow"), space = "Lab")
ColSim2 <- densCols(x = sim2$CCN4, y = sim2$Macrophages_sc_lg, nbin = 256, colramp = Lab.palette)
points(sim2$CCN4, sim2$Macrophages_sc_lg, col = ColSim2)
simData2 = data.frame(x = sim2$CCN4, y = sim2$Macrophages_sc_lg)
abline(coef(lm(y ~ x, data = simData2)), col = "red", lwd = 3)
abline(coef(lm(y ~ x, data = simData2)), col = "orange", lwd = 2)
abline(coef(lm(y ~ x, data = simData)), col = "blue", lwd = 2)

minX <- min(test2$CCN4)
maxX <- max(test2$CCN4)
minY <- min(test2$Macrophages_sc_lg)
maxY <- max(test2$Macrophages_sc_lg)

points( (test2$CCN4[test2$Cancer == 1] - minX)/(maxX - minX), (test2$Macrophages_sc_lg[test2$Cancer == 1] - minY)/(maxY - minY), col = "black")
points( (test2$CCN4[test2$Cancer == 0] - minX)/(maxX - minX), (test2$Macrophages_sc_lg[test2$Cancer == 0] - minY)/(maxY - minY), pch = 19, col = "black")

dev.off()

