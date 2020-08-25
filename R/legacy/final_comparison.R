library(tidyverse)
library(smacof)
library(psych)
library(Rtsne)
countryData <- read_csv("onlycountry_data_final.csv")

mat <- as.matrix(countryData[,c(4,6,7,12,16,19,23)])
corMat <- cor(mat)
pca <- FactoMineR::PCA(mat, scale. = F)
Mydissimilarity_10values_c <- sim2diss(corMat, method = "corr", to.dist = T)
MDS_ordinal_10values_c <- mds(Mydissimilarity_10values_c, type="ordinal")
MDSmod_10values_c <- mds(Mydissimilarity_10values_c, ndim = 2, type = c("ordinal"))
plot(MDSmod_10values_c, main = "MDS", xlab = "Dimension 1", ylab = "Dimension 2")


mat2 <- mat[,c(4,6,7,12,16,19,23)]

tsne <- Rtsne(mat2, dims = 2, pca=TRUE, perplexity = max(1,floor(nrow(mat2)/3)-1), theta=0, verbose=TRUE, max_iter = 500, check_duplicates = F)
tsne3D <- Rtsne(mat2, dims = 3, pca=TRUE, perplexity = max(1,floor(nrow(mat2)/3)-1), theta=0, verbose=TRUE, max_iter = 500, check_duplicates = F)

d_tsne_1 <- as.data.frame(tsne$Y) 
d_tsne_1 <- cbind(d_tsne_1, rowname = countryData$country)
d_tsne_1_original <-  d_tsne_1
fit_cluster_kmeans <- fpc::kmeansruns(scale(d_tsne_1[,-3]),krange=2:7,critout=F,runs=5,criterion="ch")
colpal <- randomcoloR::distinctColorPalette(fit_cluster_kmeans$bestk)
d_tsne_1_original$cl_kmeans <- factor(fit_cluster_kmeans$cluster)
plot_cex = 1
plot_cex_clus = 0.3
plot_cex_main = 0.3
plot_cex_txt = 0.3
plot_cex_txt_clus = 0.2
plot_cex_lab = 0.3
plot_cex_axis = 0.3
plot(d_tsne_1_original$V1, d_tsne_1_original$V2, col=colpal[fit_cluster_kmeans$cluster], cex = plot_cex, cex.lab = plot_cex_lab, cex.axis = plot_cex_axis, cex.main = plot_cex_main, pch=20)
text(x=d_tsne_1_original$V1, y=d_tsne_1_original$V2, cex=plot_cex_txt, pos=4, labels=(d_tsne_1_original$rowname))

data3D <- as.data.frame(tsne3D$Y)
fit_cluster_kmeans3D <- fpc::kmeansruns(scale(data3D),krange=2:7,critout=F,runs=5,criterion="ch")
colpal <- randomcoloR::distinctColorPalette(fit_cluster_kmeans3D$bestk)
fig <- plot_ly(data3D, x = ~V1, y = ~V2, z = ~V3, color = as.factor(fit_cluster_kmeans3D$cluster), 
               colors = colpal,
               text = ~paste('Term:', rownames(mat2)))
fig <- fig %>% add_markers()
fig

loadings3 <- psych::fa(corMat,nfactors = 3,fm="minres",rotate = "varimax")
