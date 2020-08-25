library(igraph)
library(graphlayouts) 
library(ggraph)
library(threejs)
library("readr")
library(phangorn)
library(phylocanvas)
library(viridis)
library(e1071)
library(stringdist)
library(seqinr)

#1st layer
links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/NZ/codon/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/NZ/codon/nodes.csv"))
colnames(links) <- c("id1","id2","label")
colnames(nodes) <- c("id","token","ord")
# ds <- readRDS("global_aligned.rda")
# names<-ds$nams
nodes$genenames <- paste(nodes$id," - ",names,sep="")
tokenCount <- plyr::count(links$label)
tokenCount$x <- as.character(tokenCount$x)
tokensOut <- tokenCount[which(tokenCount$freq>(25*(nrow(nodes)/100)) | tokenCount$freq<=(1*(nrow(nodes)/100))),1]
agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
#links<-links[which(grepl("\\+1",links$label)),]
#agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")

agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})

agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1],"-L1")
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2],"-L1")
# agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
# agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
# 
colnames(agg_links) <- c("Source","Target","Label","Weight")
#write_csv(agg_links[,-3],"multiweighted.csv")

tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L1")

gg <- graph_from_data_frame(as.data.frame(agg_links),directed = T, vertices = tmpNodes[,1])
V(gg)$label <- V(gg)$name
E(gg)$weight <- agg_links$Weight

write.csv(as.data.frame(as.matrix(as_adjacency_matrix(gg,attr = "Weight"))),"gg.csv",col.names = T,row.names = T)

#2nd layer
links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/nz/meta/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/nz/meta/nodes.csv"))
colnames(links) <- c("id1","id2","label")
nodes$genenames <- paste(nodes$node_id," - ",names,sep="")
# tokenCount <- plyr::count(links$label)
# tokenCount$x <- as.character(tokenCount$x)
# tokensOut <- tokenCount[which(tokenCount$freq<(nrow(nodes)-(20*(nrow(nodes)/100)))),1]
# agg_links <- aggregate(label ~ .,links[which(links$label %in% tokensOut),],paste, collapse = ", ")
agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")
agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})

agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
# agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
# agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
# 
colnames(agg_links) <- c("Source","Target","Label","Weight")
#write_csv(agg_links[,-3],"metaweighted.csv")

tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L2")

cluslinks <- links[which(grepl("Cluster: ",links$label)),]
cluslinks$id1 <- paste0(cluslinks$id1," - ",names[cluslinks$id1],"-L2")
cluslinks$id2 <- paste0(cluslinks$id2," - ",names[cluslinks$id2],"-L2")
colnames(cluslinks) <- c("Source","Target","Label")

write_csv(cluslinks,"metaweighted_clus.csv")

g <- graph_from_data_frame(as.data.frame(cluslinks),directed = T,vertices = tmpNodes[,5])
V(g)$label <- V(g)$name
write.csv(as.data.frame(as.matrix(as_adjacency_matrix(g))),"g.csv",col.names = T,row.names = T)
#E(g)$weight <- agg_links$Weight

tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L3")

dislinks <- links[which(grepl("District: ",links$label)),]
dislinks$id1 <- paste0(dislinks$id1," - ",names[dislinks$id1],"-L3")
dislinks$id2 <- paste0(dislinks$id2," - ",names[dislinks$id2],"-L3")
colnames(dislinks) <- c("Source","Target","Label")
write_csv(dislinks,"metaweighted_district.csv")

g2 <- graph_from_data_frame(as.data.frame(dislinks),directed = T,vertices = tmpNodes[,5])
V(g2)$label <- V(g2)$name
write.csv(as.data.frame(as.matrix(as_adjacency_matrix(g2))),"g2.csv",col.names = T,row.names = T)
#E(g2)$weight <- agg_links$Weight
#


tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L4")

dislinks <- links[which(grepl("LastCountry: ",links$label)),]
dislinks$id1 <- paste0(dislinks$id1," - ",names[dislinks$id1],"-L4")
dislinks$id2 <- paste0(dislinks$id2," - ",names[dislinks$id2],"-L4")
colnames(dislinks) <- c("Source","Target","Label")
write_csv(dislinks,"metaweighted_lastcountry.csv")

g2 <- graph_from_data_frame(as.data.frame(dislinks),directed = T,vertices = tmpNodes[,5])
V(g2)$label <- V(g2)$name
write.csv(as.data.frame(as.matrix(as_adjacency_matrix(g2))),"g3.csv",col.names = T,row.names = T)


tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L5")

dislinks <- links[which(grepl("LastCity: ",links$label)),]
dislinks$id1 <- paste0(dislinks$id1," - ",names[dislinks$id1],"-L5")
dislinks$id2 <- paste0(dislinks$id2," - ",names[dislinks$id2],"-L5")
colnames(dislinks) <- c("Source","Target","Label")
write_csv(dislinks,"metaweighted_lastcity.csv")

g2 <- graph_from_data_frame(as.data.frame(dislinks),directed = T,vertices = tmpNodes[,5])
V(g2)$label <- V(g2)$name
write.csv(as.data.frame(as.matrix(as_adjacency_matrix(g2))),"g4.csv",col.names = T,row.names = T)

tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L6")

dislinks <- links[which(grepl("OnsetDate: ",links$label)),]
dislinks$id1 <- paste0(dislinks$id1," - ",names[dislinks$id1],"-L6")
dislinks$id2 <- paste0(dislinks$id2," - ",names[dislinks$id2],"-L6")
colnames(dislinks) <- c("Source","Target","Label")
write_csv(dislinks,"metaweighted_onset.csv")

g2 <- graph_from_data_frame(as.data.frame(dislinks),directed = T,vertices = tmpNodes[,5])
V(g2)$label <- V(g2)$name
write.csv(as.data.frame(as.matrix(as_adjacency_matrix(g2))),"g5.csv",col.names = T,row.names = T)


# nucleotide layer
links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/NZ/nucleotides/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/NZ/nucleotides/nodes.csv"))
colnames(links) <- c("id1","id2","label")
nodes$genenames <- paste(nodes$id," - ",names,sep="")
tokenCount <- plyr::count(links$label)
tokenCount$x <- as.character(tokenCount$x)
tokensOut <- tokenCount[which(tokenCount$freq>(25*(nrow(nodes)/100)) | tokenCount$freq<=(1*(nrow(nodes)/100))),1]
agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
#links<-links[which(grepl("\\+1",links$label)),]
#agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")

agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})

agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1],"-L7")
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2],"-L7")
# agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
# agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
# 
colnames(agg_links) <- c("Source","Target","Label","Weight")

tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id,"-L7")
# 
#write_csv(agg_links[,-3],"nucleoweighted.csv")

gnu <- graph_from_data_frame(as.data.frame(agg_links),directed = T, vertices = tmpNodes[,1])
V(gnu)$label <- V(gnu)$name
E(gnu)$weight <- agg_links$Weight
write.csv(as.data.frame(as.matrix(as_adjacency_matrix(gnu,attr = "Weight"))),"gnu.csv",col.names = T,row.names = T)


#### the tree ------
links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/only-country_data/Taiwan/nucleotides/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/only-country_data/Taiwan/nucleotides/nodes.csv"))
colnames(links) <- c("id1","id2","label")
tmp <- readRDS(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/only-country_data/Taiwan_aligned.rda"))
nodes$genenames <- paste(nodes$id," - ",tmp$nam,sep="")
names <- tmp$nam
tokenCount <- plyr::count(links$label)
tokenCount$x <- as.character(tokenCount$x)
tokensOut <- tokenCount[which(tokenCount$freq>(100*(nrow(nodes)/100)) | tokenCount$freq<=(0*(nrow(nodes)/100))),1]
agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
#links<-links[which(grepl("\\+1",links$label)),]
#agg_links <- aggregate(label ~ .,links,paste, collapse = ", ")

agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})

agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
# agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
# agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
# 

# biased <- computeCodonBias(prep$seq)
# agg_links$weight <- agg_links$weight - unlist(lapply(agg_links$label,function(x){
#   sum(stringr::str_count(x,biased))
# }))

colnames(agg_links) <- c("Source","Target","Label","Weight")
#write_csv(agg_links[,-3],"multiweighted.csv")

tmpNodes <- nodes
tmpNodes$id <- tmpNodes$genenames
tmpNodes$id <- paste0(tmpNodes$id)

gg <- graph_from_data_frame(as.data.frame(agg_links),directed = T, vertices = tmpNodes[,1])
V(gg)$label <- V(gg)$name
E(gg)$weight <- agg_links$Weight
gp<-as.phylo(as.hclust(walktrap.community(gg)))
gp$edge.length <- computeBranchLength(gp,nodes,tokensOut,distmetric = 'set')
#saveRDS(gp$edge.length,"edge_length_global_hamming_25_1.rda")
wt <- walktrap.community(gg)

phycanv_gg <- phylocanvas(gp,treetype = "rectangular", alignlabels = F, width = 1000, height = 1000, textsize = 10, nodesize = 10, showcontextmenu = T, showscalebar = T, showhistory = F)
n <- max(wt$membership)
col_vector = viridis_pal(option = "D")(n)
for (nodename in wt$names) {
  nodename1 <- gsub("\\s","_",nodename)
  #print(nodename)
  phycanv_gg$x$nodestyles[[nodename1]]$colour <- col_vector[wt$membership[which(wt$names==nodename)]]
  #phycanv <- style_node(phycanv, nodeid = nodename, color=col_vector[wt$membership[which(wt$names==nodename)]], fillcolor=col_vector[wt$membership[which(wt$names==nodename)]])
}
phycanv_gg
write.tree(gp,"nz_tree_nucleotides_tokens_25_1_setmetric.nwk")

computeBranchLength <- function(phy,ddgnodes,tokensOut=NULL,readingframe=0,distmetric='set'){
  print("#### compute branches ####")
  # init node tokens
  nodToks <- data.frame(id=unique(c(phy$edge[,1],phy$edge[,2])),stringsAsFactors = F)
  nodToks$toks <- ""
  nodToks$order <- 0
  iter <- which(nodToks$id %in% 1:length(phy$tip.label))
  
  nodToks[iter,2] <- unlist(lapply(as.list(iter),function(x){
    leafLabel <- phy$tip.label[nodToks$id[x]]
    currentToks <- as.character(ddgnodes[which(ddgnodes$genenames==leafLabel),2])
    currentToks
  }))
  # get tree as graph and measure leaf depths from root
  ggphy <- graph.data.frame(phy$edge, directed=TRUE)
  leaves = V(ggphy)$name[which(degree(ggphy, mode = "out") == 0)]
  leaf_lengths <- all_simple_paths(ggphy,as.character(length(phy$tip.label)+1),leaves)
  initRun <- T
  while(length(leaf_lengths) > 0){
    depths <- unlist(lapply(leaf_lengths,length))
    names(depths) <- 1:length(depths)
    sorted_depths <- as.numeric(names(sort(depths,decreasing = T)))
    nextLevelLeaves <- unique(unlist(lapply(as.list(sorted_depths),function(x){
      
      current <- as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))]
      parent <- phy$edge[which(phy$edge[,2]==current),1]
      theParent <- which(nodToks$id == parent)
      theCurrent <- which(nodToks$id == current)
      
      if(initRun){
              currentToks <- str_split_fixed(as.character(ddgnodes[which(ddgnodes$genenames==phy$tip.label[as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))]]),2]),", ",n = Inf)[1,]
              currentOrder <- as.numeric(str_split_fixed(phy$tip.label[as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))]]," - ",n = Inf)[1,1])
      } else {
              currentToks <- str_split_fixed(as.character(nodToks$toks[which(nodToks$id == as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))])]),", ",n = Inf)[1,]
              currentOrder <- nodToks$order[theCurrent]
      }
      parentToks <- str_split_fixed(as.character(nodToks$toks[theParent]),", ",n = Inf)[1,]
      parentOrder <- nodToks$order[theParent]
      if(currentOrder < parentOrder | parentOrder==0){
        nodToks$order[theParent] <<- currentOrder
        combinedToks <- currentToks
      } else{
        combinedToks <- parentToks
      }
      nodToks$order[theCurrent] <<- currentOrder
      # if(distmetric=='set'){
      #   combinedToks <- base::union(currentToks,parentToks)
      # }

      parentToks <- paste(combinedToks, collapse = ", ")
      nodToks$toks[theParent] <<- parentToks
      nodToks$id[theParent]
      #
      #list(as.character(nodToks$id[which(nodToks$id == parent)]),which(nodToks$id == parent),parentToks)
    })))
    leaf_lengths <- all_simple_paths(ggphy,as.character(length(phy$tip.label)+1),as.character(nextLevelLeaves))
    initRun <- F
  }
  rm(depths)
  rm(ggphy)
  rm(sorted_depths)
  rm(nextLevelLeaves)
  gc()
  
  nodToks$toks <- as.character(nodToks$toks)
  #write_csv(as.data.frame(nodToks$order),"tmp.csv")
  
  dists <- unlist(apply(as.matrix(phy$edge),1,function(x){
    sibOne <- x[1]
    sibTwo <- x[2]
    sibOneToks <- str_split_fixed(as.character(nodToks$toks[which(nodToks$id == sibOne)]),", ",n = Inf)[1,]
    sibTwoToks <- str_split_fixed(as.character(nodToks$toks[which(nodToks$id == sibTwo)]),", ",n = Inf)[1,]

    
    if(distmetric=='set'){
      #percentage
      dist <- 1-(length(base::intersect(sibOneToks,sibTwoToks))/length(sibOneToks))
      #absolute
      #dist <- (length(base::setdiff(sibOneToks,sibTwoToks)))
    } else  if(distmetric=='order'){
      #based on order
      #if(nodToks$order[which(nodToks$id == sibOne)]==nodToks$order[which(nodToks$id == sibTwo)]){
      #  dist <- 0
      #} else{
      #  dist <- max(nodToks$order[which(nodToks$id == sibOne)],nodToks$order[which(nodToks$id == sibTwo)])
      #}
      
      dist <- abs(nodToks$order[which(nodToks$id == sibOne)] - nodToks$order[which(nodToks$id == sibTwo)])
    } else  if(distmetric=='hamming'){
      #based on hamming
      dist <- hamming.distance(sibOneToks,sibTwoToks)
      # if(dist>5000){
      #   print(base::setdiff(sibOneToks,sibTwoToks))
      # }
      #print(hamming.distance(sibOneToks,sibTwoToks))
      #print(nchar(paste(sibOneToks,collapse="")))
      #print(nchar(paste(sibTwoToks,collapse="")))
      #Biostrings::compareStrings(paste(sibOneToks,collapse=""),paste(sibTwoToks,collapse=""))
    }
    #
    #print(dist)
    
    
    
    
    if(is.nan(dist)) dist <- 0
    #dist <- 1-(length(base::union(sibOneToks,sibTwoToks))/length(sibOneToks))
    dist
  }))
  #range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  #return(scales::rescale(dists, to=c(0,1)))
  return(dists)
}

computeCodonBias <- function(prep){
  print("#### compute bias ####")
  prep <- as.data.frame(prep)
  colnames(prep) <- c("seq")
  rscus <- matrix(ncol = 64,nrow=0)
  
  for(i in 1:nrow(prep)){
    cdn <- s2c(tolower(prep$seq[i]))
    rscus <- rbind(rscus,uco(cdn,index = "rscu"))
  }
  rm(prep)
  gc()
  
  rscus.means <- colMeans(rscus)
  rscus.high <- names(which(colMeans(rscus) > 1.6))
  rscus.low <- names(which(colMeans(rscus) < 0.5))
  
  return(rscus.high)
}

#### TESTING ------
#### 
library(multinet)
net <- ml_empty()
add_layers_ml(net, c("l1","l2","l3","l4"),c(T,T,T,T))

vertices <- data.frame(V(gg)$name,"l1")
add_vertices_ml(net,vertices)
vertices <- data.frame(V(g)$name,"l2")
add_vertices_ml(net,vertices)
vertices <- data.frame(V(g2)$name,"l3")
add_vertices_ml(net,vertices)
vertices <- data.frame(V(gnu)$name,"l4")
add_vertices_ml(net,vertices)

edges <- data.frame(a=get.edgelist(gg)[,1],b="l1",c=get.edgelist(gg)[,2],d="l1")
add_edges_ml(net,edges)
edges <- data.frame(a=get.edgelist(g)[,1],b="l2",c=get.edgelist(g)[,2],d="l2")
add_edges_ml(net,edges)
edges <- data.frame(a=get.edgelist(g2)[,1],b="l3",c=get.edgelist(g2)[,2],d="l3")
add_edges_ml(net,edges)
edges <- data.frame(a=get.edgelist(gnu)[,1],b="l4",c=get.edgelist(gnu)[,2],d="l4")
add_edges_ml(net,edges)

# inter layer edges
#edges <- data.frame(a=sort(V(gg)$name),b="l1",c=sort(V(g)$name),d="l2")
#add_edges_ml(net,edges)


c <- clique_percolation_ml(net)


plot(net, vertex.labels.cex=.4, com=c, layers=c("l1","l2","l3","l4"),vertex.shape=23)
plot(net, vertex.labels.cex=.3,edge.type=3, grid=c(2,2), layers=c("l1","l2","l3","l4"))



gnup<-as.phylo(as.hclust(walktrap.community(gnu)))
gnup$edge.length <- computeBranchLength(gnup,nodes,tokensOut)
wt <- walktrap.community(gnu)
phycanv_gnu <- phylocanvas(gnup,treetype = "rectangular", alignlabels = F, width = 1000, textsize = 10, nodesize = 10, showcontextmenu = T, showscalebar = T)
n <- max(wt$membership)
col_vector = viridis_pal(option = "D")(n)
for (nodename in wt$names) {
  nodename1 <- gsub("\\s","_",nodename)
  #print(nodename)
  phycanv_gnu$x$nodestyles[[nodename1]]$colour <- col_vector[wt$membership[which(wt$names==nodename)]]
  #phycanv <- style_node(phycanv, nodeid = nodename, color=col_vector[wt$membership[which(wt$names==nodename)]], fillcolor=col_vector[wt$membership[which(wt$names==nodename)]])
}


EL  = matrix(get.edgelist(gg), ncol=2)
EL1 = matrix(get.edgelist(g), ncol=2)
EL2 = matrix(get.edgelist(g2), ncol=2)
EL3 = matrix(get.edgelist(gnu), ncol=2)

ELU = rbind(EL, EL1, EL2, EL3)
# ELU = rbind(ELU,matrix(c(unique(c(unique(EL[,1]),unique(EL[,2]))),
#                          unique(c(unique(EL1[,1]),unique(EL1[,2]))),
#                          unique(c(unique(EL2[,1]),unique(EL2[,2]))),
#                          unique(c(unique(EL1[,1]),unique(EL1[,2]))),
#                          unique(c(unique(EL2[,1]),unique(EL2[,2]))),
#                          unique(c(unique(EL3[,1]),unique(EL3[,2])))),
#                        ncol = 2))

GU = graph_from_edgelist(ELU, directed=FALSE)

V(GU)$lvl <- 1
V(GU)$lvl[which(grepl("-L2",V(GU)$name))] <- 2
V(GU)$lvl[which(grepl("-L3",V(GU)$name))] <- 3
V(GU)$lvl[which(grepl("-L4",V(GU)$name))] <- 4


#xy <- layout_as_multilevel(GU,type = "all", alpha = 25, beta = 45)

xyz <- layout_as_multilevel(GU,type = "all",
                            FUN1 = layout_as_backbone,
                            FUN2 = layout_with_stress,
                            project2D = FALSE)
GU$layout <- xyz
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
V(GU)$color <- c(sample(color, 5))[V(GU)$lvl]
V(GU)$vertex.label <- V(GU)$name

graphjs(GU, bg="white", vertex.shape="sphere",vertex.label = V(GU)$name)

library(phylocanvas)
library(dplyr)
library(ape)
plot(as.phylo(as.hclust(walktrap.community(g))))
plot(as.phylo(as.hclust(walktrap.community(gg))))


phycanv_gg
phycanv_gnu




# tesnors
library(rTensor)
dims <- c(ncol(as_adjacency_matrix(g)),ncol(as_adjacency_matrix(g)),3)

arr <- array(c(c(as.matrix(as_adjacency_matrix(g,attr = "weight"))),c(as.matrix(as_adjacency_matrix(gg,attr = "weight"))),c(as.matrix(as_adjacency_matrix(gnu,attr = "weight")))), dims)

z <- as.tensor(arr)
eDecomp <- eigen(z[,,1]@data,symmetric = F)
plot(density(eDecomp$values))
for(i in 2:3){
  eDecomp <- eigen(z[,,i]@data,symmetric = F)
  if(is.complex(eDecomp$values)){
    #plot(Re(eDecomp$values),Im(eDecomp$values))
  } else{
    lines(density(eDecomp$values))
    #plot(ecdf(eDecomp$values))
  }
}

tuckerD <- tucker(z,ranks=c(90,90,1))
coreT <- as.matrix(tuckerD$Z[,,]@data)
eDecomp <- eigen(coreT,symmetric = F)
plot(density(eDecomp$values),ylim=c(0,0.00002))
for(i in 89:2){
  tuckerD <- tucker(z,ranks=c(i,i,1))
  tuckerD$conv 
  tuckerD$norm_percent
  #plot(tuckerD$all_resids)
  
  coreT <- as.matrix(tuckerD$Z[,,]@data)
  eDecomp <- eigen(coreT,symmetric = F)
  if(is.complex(eDecomp$values)){
    #plot(Re(eDecomp$values),Im(eDecomp$values))
  } else{
    lines(density(eDecomp$values))
    #plot(ecdf(eDecomp$values))
  }
}


tsvdD <- t_svd(z)
hosvdD <-hosvd(z)
hosvdD$fnorm_resid

library("mpoly")
library(pracma)
charpoly(z[,,1]@data)
