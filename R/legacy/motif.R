library(igraph)
library(tidyverse)

setwd("~/OneDrive - Victoria University of Wellington - STAFF/RScripts/2020-03-09-epidemic-dashboard/")
folders <- c("Saudi_Arabia","Austria","France","Belgium","Japan","Colombia","Russia","Taiwan","New_Zealand")
theData <- "only-country_data"

mot3_complete <- list()
mot4_complete <- list()

mot3_complete[["20"]] <- as.data.frame(matrix(0,nrow=9,ncol=16))
mot3_complete[["100"]] <- as.data.frame(matrix(0,nrow=9,ncol=16))
mot4_complete[["20"]] <-  as.data.frame(matrix(0,nrow=9,ncol=218))
mot4_complete[["100"]] <-  as.data.frame(matrix(0,nrow=9,ncol=218))
for(j in c(20,100)){
  idx <- 1
  for(i in folders){
    #1st layer
    links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/",theData,"/",i,"/codon/links.csv"))
    nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/",theData,"/",i,"/codon/nodes.csv"))
    colnames(links) <- c("id1","id2","label")
    colnames(nodes) <- c("id","token","ord")
    
    countryDs <- readRDS(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/",theData,"/",i,"_aligned.rda"))
    names <- countryDs$nam
    nodes$genenames <- paste(nodes$id," - ",names,sep="")
    tokenCount <- plyr::count(links$label)
    tokenCount$x <- as.character(tokenCount$x)
    tokensOut <- tokenCount[which(tokenCount$freq>(j*(nrow(nodes)/100)) | tokenCount$freq<=(0*(nrow(nodes)/100))),1]
    agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
    #links<-links[which(grepl("\\+1",links$label)),]
    #agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")
    
    agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})
    
    agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1],"-L0")
    agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2],"-L0")
    # agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
    # agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
    # 
    colnames(agg_links) <- c("Source","Target","Label","Weight")
    #write_csv(agg_links[,-3],"multiweighted.csv")
    
    tmpNodes <- nodes
    tmpNodes$id <- tmpNodes$genenames
    tmpNodes$id <- paste0(tmpNodes$id,"-L0")
    
    gg <- graph_from_data_frame(as.data.frame(agg_links),directed = T, vertices = tmpNodes[,1])
    V(gg)$label <- V(gg)$name
    E(gg)$weight <- agg_links$Weight
    
    mot3 <- motifs(gg,3)
    mot3[which(is.na(mot3))] <- 0
    mot3_complete[[as.character(j)]][idx,] <- mot3
    mot4 <- motifs(gg,4)
    mot4[which(is.na(mot4))] <- 0
    mot4_complete[[as.character(j)]][idx,] <- mot4
    
    count_motifs(gg,3)
    count_motifs(gg,4)
    idx <- idx + 1
  }
}


mot3Mat <- mot3_complete[["20"]]
mot3Mat <- mot3Mat[,which(colSums(mot3Mat)>0)]
mat2 <- mot3Mat

mot4Mat <- mot4_complete[["20"]]
mot4Mat <- mot4Mat[,which(colSums(mot4Mat)>0)]
mat2 <- mot4Mat

mot3Mat <- as.data.frame(mot3Mat)
mot3Mat$ident <- rownames(mot3Mat)

pivot_longer(mot3Mat,cols = -ident) %>%
ggplot(.) +
  aes(x=ident,y=value, fill = name) +
  geom_bar(stat = "identity", position = "dodge")

mot4Mat <- as.data.frame(mot4Mat)
mot4Mat$ident <- rownames(mot4Mat)

pivot_longer(mot4Mat,cols = -ident) %>%
  ggplot(.) +
  aes(x=ident,y=value, fill = name) +
  geom_bar(stat = "identity", position = "dodge")


tsne <- Rtsne(mat2, dims = 2, pca=TRUE, perplexity = max(1,floor(nrow(mat2)/3)-1), theta=0, verbose=TRUE, max_iter = 500, check_duplicates = F)

d_tsne_1 <- as.data.frame(tsne$Y) 
d_tsne_1 <- cbind(d_tsne_1, rowname = folders)
d_tsne_1_original <-  d_tsne_1
fit_cluster_kmeans <- fpc::kmeansruns(d_tsne_1[,-3],krange=2:7,critout=F,runs=5,criterion="ch")
colpal <- randomcoloR::distinctColorPalette(fit_cluster_kmeans$bestk)
d_tsne_1_original$cl_kmeans <- factor(fit_cluster_kmeans$cluster)
plot_cex = 1
plot_cex_clus = 0.3
plot_cex_main = 0.3
plot_cex_txt = 0.5
plot_cex_txt_clus = 0.2
plot_cex_lab = 0.3
plot_cex_axis = 0.3
plot(d_tsne_1_original$V1, d_tsne_1_original$V2, col=colpal[fit_cluster_kmeans$cluster], cex = plot_cex, cex.lab = plot_cex_lab, cex.axis = plot_cex_axis, cex.main = plot_cex_main, pch=20)
text(x=d_tsne_1_original$V1, y=d_tsne_1_original$V2, cex=plot_cex_txt, pos=4, labels=(d_tsne_1_original$rowname))

