install.packages("bnlearn")
install.packages("BiocManager")
install.packages("MASS")
BiocManager::install("Rgraphviz") 

library(parallel)
library(colorspace)
library(stats) 
library(bnlearn)
library(lattice) 
library(Rgraphviz)

library(MASS)
setwd("~/Documents/Publications/BayesNetwork2/R")
rm(list = ls())

# Making compute cluster
cl = makeCluster(6, type = "SOCK")
rand = clusterEvalQ(cl, runif(10))

fileName <- "BRCA_digitalcytometry_tot_noise"
test <- read.csv(paste("./", fileName, ".csv", sep = ""), head = TRUE)

file <- "BRCA_hybrid_CIBERSORT"

DM <- c("Cancer", "CD4Tcell_sc_lg", "Neutrophils_lg", "Endothelial.cells_lg", "CAF_lg", "T.cells.CD8_lg", "NK.cells_lg", "Macrophages_sc_lg", 
        "B.cells.naive_lg", "proliferation", "Epithelial", "Mesenchymal", 
        "CCN4", "pM0", "pM1", "pM2")

test2 <- test[, DM]
test2$Cancer <- as.numeric(test2$Cancer)

# Use interval to get seed network and interval to refine the structure and infer parameters
tmp <- discretize(test2, method = "interval", breaks = c(2,rep(15,15)))

dtest2 <- as.data.frame(lapply(tmp, function(x) (as.numeric(x) - 1)/max(as.numeric(x) -1)))
dtest2 <- dtest2[complete.cases(dtest2),]

corrcoef <- cor(sapply(dtest2, as.numeric))
colnames(corrcoef) <- colnames(dtest2)
rownames(corrcoef) <- colnames(dtest2)

# Test for how normally distributed the variables are:
pdf("BRCA-HybridVariableDist-Seed-Apr20.pdf", width = 12, height = 8)
# Test for how normally distributed the variables are:
par(mfrow = c(2, 4), mar = c(4, 2, 2, 2))
for (var in DM) 
{
  x = test[, var]
  hist(x, prob = TRUE, xlab = var, ylab = "", main = "", col = "ivory")
  lines(density(x), adj = 1, lwd = 2, col = "tomato")
  #curve(dnorm(x, mean = mean(x), sd = sd(x)), from = min(x), to = max(x), add = TRUE, lwd = 2, col = "steelblue")
}

par(mfrow = c(2, 4), mar = c(4, 2, 2, 2))
for (var in names(dtest2)) 
{
  x = as.numeric(dtest2[, var])
  truehist(x, prob = TRUE, nbins = 15, h = 1/14, xlab = paste("Discretized: ", var), ylab = "", main = "", col = "ivory")
  lines(density(x), adj = 1, lwd = 2, col = "tomato")
#  curve(dnorm(x, mean = mean(x), sd = sd(x)), from = min(x), to = max(x), add = TRUE, lwd = 2, col = "steelblue")
}
dev.off()

# Plot correlation structure
library(gplots)
pdf("BRCA-Hybric-CorrelationStructure-April20.pdf", width = 9, height = 8)
  rho = cor(sapply(dtest2, as.numeric))
  #rho = cor(dtest2)
  palette.breaks = seq(-1, 1, 0.1)
  par(oma = c(2, 2, 2, 1))
  heatmap.2(rho, scale = "none", trace = "none", col = "topo.colors", revC = TRUE, breaks = palette.breaks)
dev.off()

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

testD <- dtest2
#Use a smaller number of replicates, like R = 100, for troubleshooting.
nboot = 10000 
arstr <- boot.strength(testD, R = nboot, algorithm = "iamb", cluster = cl, algorithm.args = list(blacklist = CCN4black))#, whitelist = CCN4white))
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)

CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

#Need to include whether the particular edge included in the network has a negative or a positive correlation
write.csv(cbind(arc.strength(test.iamb, testD), CorSign), paste(file, "Bootstrap", nboot, "_iamb.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "mmpc", cluster = cl) # undirected
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)
CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(ars, CorSign), paste(file, "Bootstrap", nboot, "_mmpc.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "si.hiton.pc", cluster = cl) #undirected
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)
CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(ars, CorSign), paste(file, "Bootstrap", nboot, "_si_hiton_pc.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "aracne", cluster = cl) # undirected
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)
CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(ars, CorSign), paste(file, "Bootstrap", nboot, "_aracne.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "iamb.fdr", cluster = cl, algorithm.args = list(blacklist = CCN4black))#, whitelist = CCN4white))
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)

CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(arc.strength(test.iamb, testD), CorSign), paste(file, "Bootstrap", nboot, "_iamb_fdr.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "pc.stable", cluster = cl, algorithm.args = list(blacklist = CCN4black))#, whitelist = CCN4white))
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)

CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

HLarcs <- ars[CorSign == "-",]
BRCA_str = arc.strength(test.iamb, data = testD)

strength.plot(test.iamb, BRCA_str, shape = "ellipse", highlight = list(arcs = HLarcs))

write.csv(cbind(arc.strength(test.iamb, testD), CorSign), paste(file, "Bootstrap", nboot, "_pc_stable.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "hc", cluster = cl, algorithm.args = list(blacklist = CCN4black))#, whitelist = CCN4white))
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)

CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(arc.strength(test.iamb, testD), CorSign), paste(file, "Bootstrap", nboot, "_hc.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "tabu", cluster = cl, algorithm.args = list(blacklist = CCN4black))#, whitelist = CCN4white))
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)

CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(arc.strength(test.iamb, testD), CorSign), paste(file, "Bootstrap", nboot, "_tabu.csv", sep = ' '), row.names=FALSE)

arstr <- boot.strength(testD, R = nboot, algorithm = "mmhc", cluster = cl, algorithm.args = list(blacklist = CCN4black))#, whitelist = CCN4white))
test.iamb <- averaged.network(arstr)
ars <- arcs(test.iamb)

CorSign <- rep("+", nrow(ars))

for(b in 1:nrow(ars))
{
  CorSign[b] <- ifelse(corrcoef[match(ars[b,1], colnames(corrcoef)), match(ars[b,2], colnames(corrcoef))] >0, "+", "-")
}

write.csv(cbind(arc.strength(test.iamb, testD), CorSign), paste(file, "Bootstrap", nboot, "_mmhc.csv", sep = ' '), row.names=FALSE)

# stop the cluster when we are done
stopCluster(cl)

BN_mmpc <- read.csv(file = paste(file, "Bootstrap 10000 _mmpc.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_mmpc) <- c("from", "to", "mmpc_CorSign")
BN_aracne <- read.csv(file = paste(file, "Bootstrap 10000 _aracne.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_aracne) <- c("from", "to", "aracne_CorSign")
BN_hiton <- read.csv(file = paste(file, "Bootstrap 10000 _si_hiton_pc.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_hiton) <- c("from", "to", "hiton_CorSign")
BN_iamb <- read.csv(file = paste(file, "Bootstrap 10000 _iamb.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_iamb) <- c("from", "to", "iamb_strength", "iamb_CorSign")
BN_iambfdr <- read.csv(file = paste(file, "Bootstrap 10000 _iamb_fdr.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_iambfdr) <- c("from", "to", "iambfdr_strength", "iambfdr_CorSign")
BN_tabu <- read.csv(file = paste(file, "Bootstrap 10000 _tabu.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_tabu) <- c("from", "to", "tabu_strength", "tabu_CorSign")
BN_mmhc <- read.csv(file = paste(file, "Bootstrap 10000 _mmhc.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_mmhc) <- c("from", "to", "mmhc_strength", "mmhc_CorSign")
BN_hc <- read.csv(file = paste(file, "Bootstrap 10000 _hc.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_hc) <- c("from", "to", "hc_strength", "hc_CorSign")
BN_pcstable <- read.csv(file = paste(file, "Bootstrap 10000 _pc_stable.csv", sep = ' '), head=TRUE, sep = ",", fill = TRUE, colClasses = c("character"))
names(BN_pcstable) <- c("from", "to", "pc_stable_strength", "pcstable_CorSign")

tmp <- merge(BN_mmpc, BN_aracne, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_hiton, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_iamb, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_iambfdr, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_tabu, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_mmhc, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_hc, by = c("from", "to"), all = TRUE)
tmp <- merge(tmp, BN_pcstable, by = c("from", "to"), all = TRUE)

write.csv(tmp, paste(file, "Bootstrap", nboot, "_BN_Summary_Int.csv", sep = ' '), row.names=FALSE)
