library("Biostrings")
library(stringr)
library(lubridate)
library(tesseract)
library(readr)
setwd("~/OneDrive - Victoria University of Wellington - STAFF/RScripts/2020-03-09-epidemic-dashboard/")
ref_genome <- as.character(read_csv("ref_genome.txt",col_names = F)[,1])



fils <- list.files("NZCovid/nz_data/genomes",full.names = T)

prep <- data.frame(seq=character(0),ord=numeric(0),nam=character(0),metatoks=character(0),stringsAsFactors = F)

names <- data.frame(nam=character(0),stringsAsFactors = F)
yVec <- unlist(strsplit(substr(ref_genome,55,29797),""))
#meta <- read_csv("meta.csv",col_names = F)

for(fil in fils){
  set <- readDNAStringSet(fil)
  nam <- gsub("NZCovid/nz_data/genomes/","",gsub(".consensus.fasta","",fil))
  seq <- substr(toString(set),55,29797)
  xVec <- unlist(strsplit(seq,""))
  dist <- hamming.distance(xVec,yVec)
  #seq <- toString(set)
  if(dist<=1000){
    if(length(which(meta$Accession==nam))>0){
      if(!is.na(meta$onsetdt[which(meta$Accession==nam)]) | !is.na(meta$Date.Collected[which(meta$Accession==nam)])){
        if(!is.na(meta$onsetdt[which(meta$Accession==nam)])){
          dat <- ifelse(is.na(lubridate::mdy(meta$onsetdt[which(meta$Accession==nam)])),toString(lubridate::ymd(meta$onsetdt[which(meta$Accession==nam)])),toString(lubridate::mdy(meta$onsetdt[which(meta$Accession==nam)])))
        } else {
          dat <- toString(meta$Date.Collected[which(meta$Accession==nam)])
        }
        prep <- rbind(prep,data.frame(seq=seq,ord=dat,nam=nam,metatoks=meta$tokenised[which(meta$Accession==nam)]))
      }
    }
  }
}
prep$nam <- as.character(prep$nam)
prep$seq <- as.character(prep$seq)
prep$ord <- as.character(prep$ord)
prep$metatoks <- as.character(prep$metatoks)

prep <- prep[order(prep$ord),]
names <- prep$nam


names <- as.character(names)

#saveRDS(as.data.frame(prep$seq),"seq_nz.rda")

write_csv(as.data.frame(prep$seq),"ordered_nz_only.txt",col_names = F)
write_csv(as.data.frame(prep$metatoks),"ordered_meta_nz.txt",col_names = F)


# now run via python TIC generator
# 
# now run R TIC postprocessing

library(readr)
library("Biostrings")
library(plotly)

library(ggthemes)
library(ggpubr)
library(latex2exp)
library(ggrepel)
library(scales)
library(GGally)
library(nonlinearTseries)
options(scipen=100000)

# plot coordinates -----

TICCoordinates <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/nz/TICCoordinates.csv"))

plot_ly(x=TICCoordinates$t,
        y=TICCoordinates$specificity,
        z =TICCoordinates$diversity, 
        type = "scatter3d",mode="markers",opacity=0.5) 

# plot info theory metrics -----

TICInfoTheory <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/nz/TICInfoTheory2.csv"))

plot(TICInfoTheory$ShannonWiener,type='l', xlab="Node index", ylab="Shannon entropy of unique codon matches")
plot(TICInfoTheory$Pielou,type='l', xlab="Node index", ylab="Pielou evenness of unique codon matches")


ggplot() +
  geom_line(aes(x=c(1:nrow(TICInfoTheory)), y=TICInfoTheory$ShannonWiener))+
  #geom_line(aes(x=c(1:nrow(TICInfoTheory)), y=TICInfoTheory$Pielou)) +
  theme_minimal() + 
  xlab("Progression stage") + ylab("Shannon entropy of unique codon identifiers") +
  #xlab("Progression stage") + ylab("Pielou evenness of unique codon identifiers") +
  NULL

# visualise TIC network ----

library(igraph)
links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/nz/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/nz/nodes.csv"))
colnames(links) <- c("id1","id2","label")
colnames(nodes) <- c("id","tokens","ord")
nodes$genenames <- paste(nodes$id," - ",names,sep="")
saveRDS(nodes,"ddnodes.rda")
#voc <- c("a","c","g","t")
#nonvoc <- letters[which(!letters %in% voc)]
#links <- links[which(sum(stringr::str_detect(str_split(links$label,"#")[[3]],nonvoc)) == 0),]

agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")
agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})

codonCount <- plyr::count(links$label)
codonCount <- codonCount[order(-codonCount$freq),]

links.dist <- links
links.dist$dist <- links.dist$id2 - links.dist$id1
links.dist.dens <- plyr::count(links.dist$dist)
links.dist.dens <- links.dist.dens[order(-links.dist.dens$freq),]
plot(links.dist.dens,pch=19)
plot(links.dist.dens,log='xy',pch=19)
plot(density(links.dist$dist))
ggplot(data = links.dist.dens, aes(x=x,y=freq)) + geom_line() + 
  scale_y_continuous(trans = "log") +
  scale_x_continuous(trans = "log")
write_csv(links.dist.dens,"link_dist.csv")

agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])

g <- graph_from_data_frame(as.data.frame(agg_links),directed = T)
#V(g)$label <- output$X2
#V(g)$name <- output$X2
adj <- as_adj(g, attr = "weight")
write_csv(as.data.frame(as.matrix(adj)),"adj.csv",col_names = F)

colnames(agg_links) <- c("Source","Target","label", "weight")
write_csv(agg_links,"weighted.csv")

agg_links.weightsonly <- agg_links[,-3]
agg_links.weightsonly <- agg_links.weightsonly[order(-agg_links.weightsonly$weight),]
write_csv(agg_links.weightsonly,"link_weights.csv")

V(g)$size <- 3
V(g)$frame.color <- "white"
V(g)$color <- "orange"
E(g)$arrow.mode <- 1
#E(g)$label <- links$label

l <- layout_with_dh(g)
#l <- layout_as_star(g)

plot(g, layout=l, edge.label = NA,edge.arrow.size = 0.05)

palf <- colorRampPalette(c("gold", "dark orange")) 
netm <- as.matrix(read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Git/transcendental-information-cascades/output/2020-02-03-2019nCoV/",ticRun,"/createTIC/TICmatrix.csv")))
heatmap(netm, Rowv = NA, Colv = NA, col = palf(100), scale="none", margins=c(10,10) )

# cluster analysis ----
wtc <- cluster_walktrap(g,merges = T,modularity = T)
clus <- membership(wtc)
clus.df <- data.frame(vir=names(clus), cluster=c(clus))

clus.df <- clus.df[order(clus.df$cluster),]
write_csv(clus.df,"clusters.csv")

head(clus.df)


allFrames <- c("clusters123.csv","clusters1.csv","clusters2.csv","clusters3.csv")

#intra cluster similarity
for(fram in allFrames){
  
  clustersDf <- read_csv(fram)
  
  allClusterSims <- list()
  for(i in sort(unique(clustersDf$Cluster))){
    #get all virs in cluster
    virs <- gsub("[[:digit:]]+ - ","",as.character(clustersDf$`Sequence ID`[which(clustersDf$Cluster==i)]))
    if(length(virs)>1){
      #permutartions
      intraSim <- list()
      virs.comb <- gtools::combinations(n=length(virs),2,virs,repeats.allowed = F)
      for(j in 1:nrow(virs.comb)){
        dna1 <- toString(output[which(output$X2==virs.comb[j,1]),4])
        dna1.len <- nchar(dna1)
        dna2 <- toString(output[which(output$X2==virs.comb[j,2]),4])
        dna2.len <- nchar(dna2)
        thresh <- min(dna1.len,dna2.len)
        comp <- compareStrings(substr(dna1,1,thresh), substr(dna2,1,thresh))
        comp.split <- unlist(strsplit(comp,""))
        dnaDiff <- length(which(comp.split=="?"))
        dnaDiff <- dnaDiff + length(which(comp.split=="-"))
        dnaDiff <- dnaDiff + length(which(comp.split=="+"))
        intraSim[[length(intraSim)+1]] <- 1 - (dnaDiff/thresh)
      }
      
      allClusterSims[[length(allClusterSims)+1]] <- intraSim
    } else{
      allClusterSims[[length(allClusterSims)+1]] <- 1
    }
  }
  
  avgSim <- lapply(allClusterSims,function(x){
    c(mean(unlist(x)),sd(unlist(x)),min(unlist(x)),max(unlist(x)))
  })
  
  clSizes <- plyr::count(clustersDf$Cluster)
  
  plotData <- data.frame(lab=sort(unique(clustersDf$Cluster)),avgSim=unlist(avgSim)[seq(from=1,to=length(unlist(avgSim)),by=4)],sdSim=unlist(avgSim)[seq(from=2,to=length(unlist(avgSim)),by=4)],minSim=unlist(avgSim)[seq(from=3,to=length(unlist(avgSim)),by=4)],maxSim=unlist(avgSim)[seq(from=4,to=length(unlist(avgSim)),by=4)],clSize=clSizes$freq)
  
  p <- ggplot(data=plotData, aes(x=as.character(lab), y=avgSim)) +
    geom_bar(stat="identity", fill="#E69F00")+
    geom_text(aes(label=round(avgSim,2)), vjust=1.6, color="white", size=3.5)+
    theme_minimal() + 
    xlab("Cluster") + ylab("Avg similarity")
  ggsave(plot = p, filename = paste0(gsub("\\.csv","_stat",fram),".png"),device = png())
  
  write_csv(plotData,paste0(gsub("\\.csv","_stat",fram),".csv"))
  
  # inter cluster similarity
  
  clusIds <- sort(unique(clustersDf$Cluster))
  clusIds.combs <- gtools::combinations(n=length(clusIds),2,clusIds,repeats.allowed = F)
  
  allInterClusterSims <- list()
  
  for(i in 1:nrow(clusIds.combs)){
    nextComb1 <- clusIds.combs[i,1]
    nextComb2 <- clusIds.combs[i,2]
    
    nextComb1.virs <- gsub("[[:digit:]]+ - ","",as.character(clustersDf$`Sequence ID`[which(clustersDf$Cluster==nextComb1)]))
    nextComb2.virs <- gsub("[[:digit:]]+ - ","",as.character(clustersDf$`Sequence ID`[which(clustersDf$Cluster==nextComb2)]))
    
    interCombs <- expand.grid(nextComb1.virs,nextComb2.virs)
    interSim <- list()
    for(j in 1:nrow(interCombs)){
      dna1 <- toString(output[which(output$X2==interCombs[j,1]),4])
      dna1.len <- nchar(dna1)
      dna2 <- toString(output[which(output$X2==interCombs[j,2]),4])
      dna2.len <- nchar(dna2)
      thresh <- min(dna1.len,dna2.len)
      comp <- compareStrings(substr(dna1,1,thresh), substr(dna2,1,thresh))
      comp.split <- unlist(strsplit(comp,""))
      dnaDiff <- length(which(comp.split=="?"))
      dnaDiff <- dnaDiff + length(which(comp.split=="-"))
      dnaDiff <- dnaDiff + length(which(comp.split=="+"))
      interSim[[length(interSim)+1]] <- 1 - (dnaDiff/thresh)
    }
    allInterClusterSims[[paste0(clusIds.combs[i,],collapse="-")]] <- interSim
  }
  
  pairNames <- names(allInterClusterSims)
  
  avgSim <- lapply(allInterClusterSims,function(x){
    c(mean(unlist(x)),sd(unlist(x)),min(unlist(x)),max(unlist(x)))
  })
  
  plotData <- data.frame(lab=pairNames,avgSim=unlist(avgSim)[seq(from=1,to=length(unlist(avgSim)),by=4)],sdSim=unlist(avgSim)[seq(from=2,to=length(unlist(avgSim)),by=4)],minSim=unlist(avgSim)[seq(from=3,to=length(unlist(avgSim)),by=4)],maxSim=unlist(avgSim)[seq(from=4,to=length(unlist(avgSim)),by=4)])
  write_csv(plotData,paste0(gsub("\\.csv","_interstat",fram),".csv"))
  p <- ggplot(data=plotData, aes(x=lab, y=avgSim)) +
    geom_bar(stat="identity", fill="#E69F00")+
    geom_text(aes(label=round(avgSim,2)), vjust=1.6, color="white", size=3.5)+
    theme_minimal() +
    xlab("Inter-cluster pair") + ylab("Avg similarity")
  ggsave(plot = p, filename = paste0(gsub("\\.csv","_interstat",fram),".png"),device = png(), width = 12)
}

# RQA ----
rq0bj <- rqa(time.series = TICCoordinates$diversity,radius = 1,time.lag = 1)
plot(rq0bj)
write_csv(as.data.frame(TICCoordinates$diversity),"div_series.csv")


# generate other networks -----
# +1
filteredLinks <- links[which(grepl("\\+1",links$label)),]
colnames(filteredLinks) <- c("id1","id2","label")
agg_links2 <- aggregate(label ~ .,filteredLinks,paste, collapse = ", ")
agg_links2$weight <- sapply(agg_links2$label,function(x){stringr::str_count(x,", ")+1})
agg_links2$id1 <- paste0(rownames(output)[agg_links2$id1]," - ",output$X2[agg_links2$id1])
agg_links2$id2 <- paste0(rownames(output)[agg_links2$id2]," - ",output$X2[agg_links2$id2])
colnames(agg_links2) <- c("Source","Target","label", "weight")
write_csv(agg_links2,"weighted_+1.csv")

agg_links2.weightsonly <- agg_links2[,-3]
agg_links2.weightsonly <- agg_links2.weightsonly[order(-agg_links2.weightsonly$weight),]
write_csv(agg_links2.weightsonly,"link_weights_+1.csv")

# +2
filteredLinks <- links[which(grepl("\\+2",links$label)),]
colnames(filteredLinks) <- c("id1","id2","label")
agg_links2 <- aggregate(label ~ .,filteredLinks,paste, collapse = ", ")
agg_links2$weight <- sapply(agg_links2$label,function(x){stringr::str_count(x,", ")+1})
agg_links2$id1 <- paste0(rownames(output)[agg_links2$id1]," - ",output$X2[agg_links2$id1])
agg_links2$id2 <- paste0(rownames(output)[agg_links2$id2]," - ",output$X2[agg_links2$id2])
colnames(agg_links2) <- c("Source","Target","label", "weight")
write_csv(agg_links2,"weighted_+2.csv")

agg_links2.weightsonly <- agg_links2[,-3]
agg_links2.weightsonly <- agg_links2.weightsonly[order(-agg_links2.weightsonly$weight),]
write_csv(agg_links2.weightsonly,"link_weights_+2.csv")

# +3
filteredLinks <- links[which(grepl("\\+3",links$label)),]
colnames(filteredLinks) <- c("id1","id2","label")
agg_links2 <- aggregate(label ~ .,filteredLinks,paste, collapse = ", ")
agg_links2$weight <- sapply(agg_links2$label,function(x){stringr::str_count(x,", ")+1})
agg_links2$id1 <- paste0(rownames(output)[agg_links2$id1]," - ",output$X2[agg_links2$id1])
agg_links2$id2 <- paste0(rownames(output)[agg_links2$id2]," - ",output$X2[agg_links2$id2])
colnames(agg_links2) <- c("Source","Target","label", "weight")
write_csv(agg_links2,"weighted_+3.csv")

agg_links2.weightsonly <- agg_links2[,-3]
agg_links2.weightsonly <- agg_links2.weightsonly[order(-agg_links2.weightsonly$weight),]
write_csv(agg_links2.weightsonly,"link_weights_+3.csv")

# compare original clustering with tandom-----
library(clusteval)
library(vcd)
origClus <- read_csv("clusters123.csv")
randomClus <- read_csv("rand-clusters.csv")

origClus <- arrange(origClus,`Sequence ID`)
randomClus <- arrange(randomClus,`Sequence ID`)

clSim <- cluster_similarity(origClus$Cluster, randomClus$Cluster, similarity = c("jaccard"))
coMem <- comembership_table(origClus$Cluster, randomClus$Cluster)

coMat <- matrix(c(coMem$n_11, coMem$n_10, coMem$n_01, coMem$n_00), nrow = 2, byrow = T)

ctable <- as.table(coMat)
fourfoldplot(ctable,
             conf.level = 0, margin = 1, main = "Confusion Matrix")

library(irr)
kappa2(cbind(as.factor(origClus$Cluster),as.factor(randomClus$Cluster)), "unweighted")

