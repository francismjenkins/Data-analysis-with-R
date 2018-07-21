# Author: Frank M. Jenkins
# Project: Data Analysis, BINF 702
# Spring 2018, Professor Solka


library(rcellminer); library(factoextra); data(molData)

nciExp <- getALLFeatureData(rcellminerData::molData)[["xai"]]

nciExp2 <- nciExp[1:50,1:18]


nciScale <- scale(na.omit(nciExp2), center = T, scale = T)

nci.pca <- princomp(na.omit(nciScale), cor = T, scores = T)

fviz_screeplot(nci.pca)
summary(nci.pca)
fviz_pca_ind(nci.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

nciScale2 <- t(nciScale)

nci.clust <- hclust(dist(nciScale2, method = "euclidean", diag = F, upper = F, p = 2), method = "average", members = NULL)

plot(nci.clust)

nciBrCNS <- nciExp[,1:10]

nciCors <- cor(nciBrCNS, use = "complete.obs", method = "spearman")
head(nciCors)

library(corrplot)

nciBrCNS <- t(nciBrCNS)
nciBrCNS <- nciBrCNS[,1:20]
nciCor2 <- cor(nciBrCNS, use = "complete.obs", method = "spearman")
corrplot(nciCor2)

nci60Mut <- getAllFeatureData(rcellminerData::molData)[["mut"]]

DNArepair <- c("SMUG1", "MYH", "LIG3", "XRCC1", "MGMT", "MSH2", "PMS1", "MLH1", "APC", "APLF", "ATM",
"CLSPN", "ERCC6", "FANCI", "FANCM", "GEN1", "HLTF", "MLH1", "POLD1", "POLE", "POLG", "POLQ", "RAD54L", 
"REV3L", "RMI1", "SLX4", "SMARCA4", "SMC2", "TP53", "WRN")
DNArepair <- intersect(DNArepair, rownames(nci60Mut))

DNAmut <- nci60Mut[DNArepair, ]

numMut <- apply(DNAmut, MARGIN = 1, FUN = sum)
sort(numMut, decreasing = TRUE)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)


mutPlot <- DNAmut[order(numMut), ]
heatmap.2(mutPlot, main = "Deleterious Mutations", col = my_palette)
