library(igraph)
library(ape)
# load our tree

links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/New_Zealand2/nucleotides/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/New_Zealand2/nucleotides/nodes.csv"))
colnames(links) <- c("id1","id2","label")
colnames(nodes) <- c("id","token","ord")
nodes$genenames <- paste(nodes$id," - ",names,sep="")
tokenCount <- plyr::count(links$label)
tokenCount$x <- as.character(tokenCount$x)
tokensOut <- tokenCount[which(tokenCount$freq>(70*(nrow(nodes)/100)) | tokenCount$freq<=(0*(nrow(nodes)/100))),1]
agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
#links<-links[which(grepl("\\+1",links$label)),]
#agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")

agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})

agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
colnames(agg_links) <- c("source","target","label","weight")
nodes$id <- nodes$genenames
g <- graph_from_data_frame(as.data.frame(agg_links),directed = T,vertices = nodes)
wtc <- cluster_walktrap(g,merges = T,modularity = T)
phy <- as.phylo(as.hclust(wtc))

#load other tree

nxtree <- read.tree("nextstrain_groups_nz-covid19-private_2020-04-29_discrete-phylogeo_278-time_tree.nwk")

phy$tip.label <- gsub("([[:digit:]]+) - hCoV-19/|([[:digit:]]+) - ","",phy$tip.label)
nxtree$tip.label <- gsub("/[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]","",gsub("NewZealand/","",unlist(lapply(strsplit(nxtree$tip.label,"\\|"),function(x){x[1]}))))
# compare trees

#treeComp <- comparePhylo(phy,nxtree, plot = T)
phy$tip.label[which(phy$tip.label %in% nxtree$tip.label)]

# library(distory)
# library(phytools)
# library(kdetrees)
# phy$edge.length<-1
# nxtree$edge.length<-1
# 
# plot(nxtree, font = 1, cex = 0.25)

#phylo.diff(phy,nxtree)
# library(dendextend)
# tanglegram(as.dendrogram(as.hclust(force.ultrametric(phy,method = "nnls"))), as.dendrogram(as.hclust(force.ultrametric(nxtree,method = "nnls"))))
# cor_bakers_gamma(phy, nxtree)


# library(graphkernels)
# K <- CalculateWLKernel(list(as.igraph(phy),as.igraph(nxtree)), 5)

# using shortest paths over the common nodes
g1 <- as.igraph(phy,directed = F)
g2 <- as.igraph(nxtree, directed = F)
commonLabels <- phy$tip.label[which(phy$tip.label %in% nxtree$tip.label)]
commonLabels.combinations <- as.data.frame(gtools::combinations(length(commonLabels),2,commonLabels))
commonLabels.combinations$distG1 <- 0
commonLabels.combinations$distG2 <- 0

commonLabels.combinations$distG1<- unlist(lapply(1:nrow(commonLabels.combinations),function(x,y){
  nextPath <- shortest_paths(g1,commonLabels.combinations[x,1],commonLabels.combinations[x,2])
  length(nextPath$vpath[[1]]$name)
},y=commonLabels.combinations))

commonLabels.combinations$distG2<- unlist(lapply(1:nrow(commonLabels.combinations),function(x,y){
  nextPath <- shortest_paths(g2,commonLabels.combinations[x,1],commonLabels.combinations[x,2])
  length(nextPath$vpath[[1]]$name)
},y=commonLabels.combinations))


hist(abs(commonLabels.combinations$distG1 - commonLabels.combinations$distG2))
hist(commonLabels.combinations$distG1 - commonLabels.combinations$distG2)
dev.off()
plot(density(commonLabels.combinations$distG1 - commonLabels.combinations$distG2))

g <- commonLabels.combinations$distG1 - commonLabels.combinations$distG2
m <- mean(g)
std <- sd(g)
hist(g , main="Path difference between common genomes in trees",xlab="x",col="green",label=TRUE,plot = TRUE, freq = F, ylim = c(0,0.15))  
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")
lines(density(g))
abline(v = m, lty = 2)
abline(v = m+std, lty = 2)
abline(v = m-std, lty = 2)
abline(v = m+2*std, lty = 2)
abline(v = m-2*std, lty = 2)

diameter(g1)
diameter(g2)

mean_distance(g1)
mean_distance(g2)








#### find optimal filter --------
out <- data.frame(i=numeric(0),j=numeric(0),m1=numeric(0),std=numeric,m2=numeric(0))
for(i in 100:20){
  for(j in 0:19){
    tokensOut <- tokenCount[which(tokenCount$freq>(i*(nrow(nodes)/100)) | tokenCount$freq<=(j*(nrow(nodes)/100))),1]
    agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
    
    agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})
    
    agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
    agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
    colnames(agg_links) <- c("source","target","label","weight")
    nodes$id <- nodes$genenames
    g <- graph_from_data_frame(as.data.frame(agg_links),directed = T,vertices = nodes)
    wtc <- cluster_walktrap(g,merges = T,modularity = T)
    phy <- as.phylo(as.hclust(wtc))
    nxtree <- read.tree("nextstrain_groups_nz-covid19-private_2020-04-29_discrete-phylogeo_278-time_tree.nwk")
    phy$tip.label <- gsub("([[:digit:]]+) - hCoV-19/|([[:digit:]]+) - ","",phy$tip.label)
    nxtree$tip.label <- gsub("/[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]","",gsub("NewZealand/","",unlist(lapply(strsplit(nxtree$tip.label,"\\|"),function(x){x[1]}))))
    phy$tip.label[which(phy$tip.label %in% nxtree$tip.label)]
    g1 <- as.igraph(phy,directed = F)
    g2 <- as.igraph(nxtree, directed = F)
    commonLabels <- phy$tip.label[which(phy$tip.label %in% nxtree$tip.label)]
    commonLabels.combinations <- as.data.frame(gtools::combinations(length(commonLabels),2,commonLabels))
    commonLabels.combinations$distG1 <- 0
    commonLabels.combinations$distG2 <- 0
    commonLabels.combinations$distG1<- unlist(lapply(1:nrow(commonLabels.combinations),function(x,y){
      nextPath <- shortest_paths(g1,commonLabels.combinations[x,1],commonLabels.combinations[x,2])
      length(nextPath$vpath[[1]]$name)
    },y=commonLabels.combinations))
    
    commonLabels.combinations$distG2<- unlist(lapply(1:nrow(commonLabels.combinations),function(x,y){
      nextPath <- shortest_paths(g2,commonLabels.combinations[x,1],commonLabels.combinations[x,2])
      length(nextPath$vpath[[1]]$name)
    },y=commonLabels.combinations))
    
    g <- commonLabels.combinations$distG1 - commonLabels.combinations$distG2
    m <- mean(g)
    std <- sd(g)
    med <- median(g)
    out <- rbind(out,data.frame(i=i,j=j,m1=m,std=std,m2=med))
    print(i)
    print(j)
    print(m)
    print(std)
    print(med)
  }
}

