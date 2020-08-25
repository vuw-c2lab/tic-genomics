library(tidyverse)
library("Biostrings")
library(stringr)
library(lubridate)
library(tesseract)
library(readr)
library(e1071)
library(stringdist)
setwd("~/OneDrive - Victoria University of Wellington - STAFF/RScripts/2020-03-09-epidemic-dashboard/")
ref_genome <- as.character(read_csv("ref_genome.txt",col_names = F)[,1])

#subset.ds <- readRDS("global_aligned_subset.rda")
prep <- readRDS("plain_global_aligned.rda")

fils <- list.files("NZCovid/nz_data/wavetwo",full.names = T)

prepNZ <- prep[which(grepl("NewZealand",prep$nams) & grepl("20VR",prep$nams)),]
nzNames <- unlist(lapply(strsplit(prepNZ$nams,"/"),function(x){
  x[2]
}))
prepNZ$nams <- nzNames

filsQuarantine <- gsub(".consensus.fasta","",list.files("NZCovid/nz_data/newgenomes",full.names = F))
filsQuarantine <- filsQuarantine[-length(filsQuarantine)]
filsQuarantine <- filsQuarantine[which(!filsQuarantine %in% nzNames)]

prepRest <- prep[which(!grepl("NewZealand",prep$nams)),]

prep <- prepNZ

whereFrom <- c("Brazil","England","USA","Australia","Belgium","Spain","China","Italy","Netherlands","Poland","India","HongKong","UnitedArabEmirates","Norway","Taiwan","Pakistan","Austria","Sweden","Scotland","Wales","France")
for(country in whereFrom){
  prepCountry <- prepRest[which(grepl(country,prepRest$nams)),]
  subsample <- unique(sort(c(sample(1:(nrow(prepCountry)),min(c(40,nrow(prepCountry))),replace = F))))
  prep <- rbind(prep,prepCountry[subsample,])
}

prep <- prep[order(prep$dats),]

names <- data.frame(nam=character(0),stringsAsFactors = F)

meta <- read_csv("NZCovid/nz_data/newmeta.csv",col_names = F)
meta$X8 <- as.character(lubridate::dmy(meta$X8))
yVec <- unlist(strsplit(substr(ref_genome,55,29797),""))
prep <- rbind(prep,data.frame(seq=paste(yVec,collapse = ""),dats="2019-12-01",nams="China/Wuhan-Hu-1"))

for(fil in filsQuarantine){
  set <- readDNAStringSet(paste0("NZCovid/nz_data/newgenomes/",fil,".consensus.fasta"))
  nam <- fil
  seq <- substr(toString(set),55,29797)

  if(nchar(seq) == length(yVec) & length(which(meta$X1==nam))>0){
    dat <- meta$X8[which(meta$X1==nam)]
    if(grepl("2020",dat) & nchar(dat)==10 & !grepl(nam,prep$nams)){
      xVec <- unlist(strsplit(seq,""))
      dist <- hamming.distance(xVec,yVec)
      #print(dist)
      #print(length(base::intersect(xVec,yVec)))/length(xVec)
      #seq <- toString(set)
      match <- regexpr("(N)\\1+",seq)
      print(attributes(match)$match.length)
      #if(dist<=5000 & max(attributes(match)$match.length < 100)){
      prep <- rbind(prep,data.frame(seq=seq,dats=dat,nams=nam))
      #}
    }
  }
}

prep <- prep[order(prep$dats),]
prep <- prep[,-2]

# for(fil in fils){
#   set <- readDNAStringSet(fil)
#   nam <- gsub("NZCovid/nz_data/wavetwo/","",gsub(".consensus.fasta","",gsub(".fasta","",gsub(".consensus","",fil))))
#   seq <- substr(toString(set),55,29797)
# 
#   if(nchar(seq) == length(yVec)){
# 
#     xVec <- unlist(strsplit(seq,""))
#     dist <- hamming.distance(xVec,yVec)
#     #print(dist)
#     #print(length(base::intersect(xVec,yVec)))/length(xVec)
#     #seq <- toString(set)
#     match <- regexpr("(N)\\1+",seq)
#     print(attributes(match)$match.length)
#     #if(dist<=5000 & max(attributes(match)$match.length < 100)){
#       prep <- rbind(prep,data.frame(seq=seq,nams=nam))
#     #}
#   }
# }

prep$nams <- as.character(prep$nams)
prep$seq <- as.character(prep$seq)

#prep <- unique(data.table::setDT(prep), by = c("nams"))

names <- prep$nams


names <- as.character(names)

#saveRDS(as.data.frame(prep$seq),"New_Zealand-seq.rda")
#saveRDS(prep,"New_Zealand_aligned.rda")

write_csv(as.data.frame(prep$seq),"plain_analysis.txt",col_names = F)
setwd("/Users/MLR/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/plain_analysis/")
saveRDS(prep,"plain_analysis_aligned.rda")


### the tree ------
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

links <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/plain_analysis/links.csv"))
nodes <- read_csv(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/plain_analysis/nodes.csv"))
colnames(links) <- c("id1","id2","label")
tmp <- readRDS(paste0("~/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/plain_analysis/plain_analysis_aligned.rda"))
nodes$genenames <- paste(nodes$id," - ",tmp$nam,sep="")
names <- tmp$nam
tokenCount <- plyr::count(links$label)
tokenCount$x <- as.character(tokenCount$x)
orderedTokenCount <- tokenCount[order(-tokenCount$freq),2]
#tokensOut <- tokenCount[which(tokenCount$freq>(100*(nrow(nodes)/100)) | tokenCount$freq<=(0*(nrow(nodes)/100))),1]
tokensOut <- tokenCount[which(tokenCount$freq>=(orderedTokenCount[which(diff(orderedTokenCount)==min(diff(orderedTokenCount)))]) | tokenCount$freq<=1),1]
#tokensOut <- tokenCount[which(tokenCount$freq>=100 | tokenCount$freq<=0),1]
agg_links <- aggregate(label ~ .,links[which(!links$label %in% tokensOut),],paste, collapse = ", ")
orderedLinkCount <- plyr::count(agg_links$id2)
orderedLinkCount <- orderedLinkCount[order(-orderedLinkCount$freq),]
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

gg <- graph_from_data_frame(as.data.frame(agg_links),directed = F, vertices = tmpNodes[,1])
V(gg)$label <- V(gg)$name
E(gg)$weight <- agg_links$Weight
wt <- cluster_fast_greedy(gg)
gp<-as.phylo(as.hclust(wt))
#gp$edge.length <- computeBranchLength(gp,nodes,tokensOut,distmetric = 'hamming',prep=prep)
#saveRDS(gp$edge.length,"edge_length_global_hamming_25_1.rda")
gp2 <- root(gp, out=1)
gp2 <- ladderize(gp2)
phycanv_gg <- phylocanvas(gp2,treetype = "rectangular", alignlabels = F, width = "100%", height = 1000, textsize = 10, nodesize = 10, showcontextmenu = T, showscalebar = T, showhistory = F)
n <- max(wt$membership)
col_vector = viridis_pal(option = "D")(n)
for (nodename in wt$names) {
  nodename1 <- gsub("\\s","_",nodename)
  #print(nodename)
  phycanv_gg$x$nodestyles[[nodename1]]$colour <- col_vector[wt$membership[which(wt$names==nodename)]]
  #phycanv <- style_node(phycanv, nodeid = nodename, color=col_vector[wt$membership[which(wt$names==nodename)]], fillcolor=col_vector[wt$membership[which(wt$names==nodename)]])
}
phycanv_gg
write.tree(gp,"plain_analysis_tree.nwk")


computeBranchLength <- function(phy,ddgnodes,tokensOut=NULL,readingframe=0,distmetric='set',prep){
  print("#### compute branches ####")
  # init node tokens
  nodToks <- data.frame(id=unique(c(phy$edge[,1],phy$edge[,2])),stringsAsFactors = F)
  nodToks$toks <- ""
  nodToks$order <- 0
  nodToks$seqs <- ""
  iter <- which(nodToks$id %in% 1:length(phy$tip.label))
  
  nodToks[iter,2] <- unlist(lapply(as.list(iter),function(x){
    leafLabel <- phy$tip.label[nodToks$id[x]]
    currentToks <- as.character(ddgnodes[which(ddgnodes$genenames==leafLabel),2])
    currentToks
  }))
  nodToks[iter,4] <- unlist(lapply(as.list(iter),function(x){
    leafLabel <- phy$tip.label[nodToks$id[x]]
    currentToks <- as.character(prep[which(ddgnodes$genenames==leafLabel),1])
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
        currentSeqs <- str_split_fixed(as.character(prep[which(ddgnodes$genenames==phy$tip.label[as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))]]),1]),"",n = Inf)[1,]
      } else {
        currentToks <- str_split_fixed(as.character(nodToks$toks[which(nodToks$id == as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))])]),", ",n = Inf)[1,]
        currentOrder <- nodToks$order[theCurrent]
        currentSeqs <- str_split_fixed(as.character(nodToks$seqs[which(nodToks$id == as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))])]),"",n = Inf)[1,]
      }
      parentToks <- str_split_fixed(as.character(nodToks$toks[theParent]),", ",n = Inf)[1,]
      parentOrder <- nodToks$order[theParent]
      parentSeqs <- str_split_fixed(as.character(nodToks$seqs[theParent]),"",n = Inf)[1,]
      if(currentOrder < parentOrder | parentOrder==0){
        nodToks$order[theParent] <<- currentOrder
        combinedToks <- currentToks
        combinedSeqs <- currentSeqs
      } else{
        combinedToks <- parentToks
        combinedSeqs <- parentSeqs
      }
      nodToks$order[theCurrent] <<- currentOrder
      # if(distmetric=='set'){
      #   combinedToks <- base::union(currentToks,parentToks)
      # }
      
      parentToks <- paste(combinedToks, collapse = ", ")
      parentSeqs <- paste(combinedSeqs, collapse = "")
      nodToks$toks[theParent] <<- parentToks
      nodToks$seqs[theParent] <<- parentSeqs
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
    } else if(distmetric=='gene'){
      aDNA <- unlist(strsplit(as.character(nodToks$seqs[which(nodToks$id == sibOne)]),""))
      bDNA <- unlist(strsplit(as.character(nodToks$seqs[which(nodToks$id == sibTwo)]),""))
      
      DNAMat <- matrix(c(aDNA,bDNA),nrow = 2,byrow = T)
      dist <- dist.dna(as.DNAbin(list(aDNA,bDNA)),model="TN93")[1]
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



#### the network-----
library(visNetwork)
dd_filt <- agg_links
gg <- graph.data.frame(dd_filt, directed=F)

# build tree from clusters
wt_b <- cluster_fast_greedy(gg)
n <- max(wt_b$membership)
gg <- graph.data.frame(dd_filt, directed=T)
col_vector = viridis_pal(option = "D")(n)
visNLinks <- dd_filt[,-3]
visNNodes <- as.data.frame(unique(c(dd_filt$Source,dd_filt$Target)))
colnames(visNNodes) <- c("id")
visNNodes$id <- as.character(visNNodes$id)
visNNodes$group <- wt_b$membership
visNNodes$label <- visNNodes$id
visNNodes$title <- paste0("<p><b>", visNNodes$id,"</b></p>")
visNNodes$value <- 100
visNNodes$size <- 100
visNNodes$color <- col_vector[wt_b$membership]
colnames(visNLinks) <- c("from","to","weight")
visNLinks$value <- visNLinks$weight
visNLinks$label <- visNLinks$weight
visNLinks$font.size <- 40
visNetwork(nodes=visNNodes ,edges = visNLinks, width="100%", height="800px") %>% 
  visEdges(arrows = "to") %>% 
  visNodes(font = list(size = 72)) %>%
  visOptions(manipulation = F, selectedBy = "group",highlightNearest = TRUE) %>% 
  visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -500)) %>%
  visInteraction(keyboard = TRUE, tooltipDelay = 0, hideEdgesOnDrag = TRUE)

