library(phylocanvas)
library(dplyr)
library(ape)
library(viridis)
phycanv <- phylocanvas(read.tree(text="(A:5,B:0.2,(C:0.3,D:0.4)E:0.5)F;"), treetype = "rectangular", alignlabels = F)


phy <- as.phylo(as.hclust(wt))

# compute edge lengths
phy$edge.length

phy$edge.length <- sample(1,length(phy$edge.length),replace = T)

tmpNodes <- strsplit(phy$tip.label[phy$edge[which(phy$edge[,2] %in% 1:21),2]], " - ")
tmpNodes <- as.numeric(unlist(lapply(tmpNodes, function(x){
  x[1]
})))

phy$edge.length[which(phy$edge[,2] %in% 1:21)] <- tmpNodes

phycanv <- phylocanvas(phy, treetype = "rectangular", alignlabels = F)

library("phangorn")
library(Biostrings)
library("ape")
library(msa)
library(phylocanvas)
fsta <- read.FASTA("gisaid_cov2020_sequences_belgium.fasta", type = "DNA")
alin <- msa::msaMuscle(prep$seq[1:2], type = "dna")
primates <- read.phyDat("exdna.txt", format="phylip", type="DNA")
dm = dist.dna(fsta)

pairwiseAlignment(DNAString(prep$seq[1]),DNAString(prep$seq[2]))

library(ape)
nextstrain.oceania <- read.tree("~/Downloads/nextstrain_ncov_oceania_tree.nwk")
plot(nextstrain.oceania)



#### filter out commons
computeBranchLength <<- function(phy,ddgnodes,tokensOut=NULL,readingframe=0){
  print("#### compute branches ####")
  # init node tokens
  nodToks <- data.frame(id=unique(c(phy$edge[,1],phy$edge[,2])))
  nodToks$toks <- ""
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
      parent <- as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))-1]
      currentToks <- ifelse (initRun,
                             str_split_fixed(as.character(ddgnodes[which(ddgnodes$genenames==phy$tip.label[as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))]]),2]),", ",n = Inf)[1,],
                             str_split_fixed(as.character(nodToks$toks[which(nodToks$id == as.numeric(leaf_lengths[[x]]$name)[length(as.numeric(leaf_lengths[[x]]$name))])]),", ",n = Inf)[1,]
      )
      theParent <- which(nodToks$id == parent)
      parentToks <- str_split_fixed(as.character(nodToks$toks[theParent]),", ",n = Inf)[1,]
      combinedToks <- base::union(currentToks,parentToks)
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
  dists <- unlist(apply(as.matrix(phy$edge),1,function(x){
    sibOne <- x[1]
    sibTwo <- x[2]
    sibOneToks <- str_split_fixed(as.character(nodToks$toks[which(nodToks$id == sibOne)]),", ",n = Inf)[1,]
    sibTwoToks <- str_split_fixed(as.character(nodToks$toks[which(nodToks$id == sibTwo)]),", ",n = Inf)[1,]
    
    #percentage
    #dist <- (length(base::setdiff(sibOneToks,sibTwoToks))/length(sibOneToks))*100
    #absolute
    dist <- (length(base::setdiff(sibOneToks,sibTwoToks)))
    if(is.nan(dist)) dist <- 0
    
    #dist <- 1-(length(base::union(sibOneToks,sibTwoToks))/length(sibOneToks))
    dist
  }))
  #range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  return(scales::rescale(dists, to=c(0,1)))
}

tokenCount <- plyr::count(links$label)
tokenCount$x <- as.character(tokenCount$x)
#tokenCount <- tokenCount[order(-tokenCount$freq),]
#tokenCount[which(tokenCount$freq<(nrow(nodes)-9)),1]

tokensOut <- tokenCount[which(tokenCount$freq<(nrow(nodes)-9)),1]

agg_links <- aggregate(label ~ .,links[which(links$label %in% tokensOut),],paste, collapse = ", ")
agg_links$weight <- sapply(agg_links$label,function(x){stringr::str_count(x,", ")+1})
agg_links$id1 <- paste0(agg_links$id1," - ",names[agg_links$id1])
agg_links$id2 <- paste0(agg_links$id2," - ",names[agg_links$id2])
colnames(agg_links) <- c("Source","Target","label", "weight")
#write_csv(agg_links,"weighted.csv")

g <- graph_from_data_frame(as.data.frame(agg_links),directed = T)
wt <- cluster_walktrap(g,merges = T,modularity = T)
phy <- as.phylo(as.hclust(wt))

phy$edge.length <- computeBranchLength(phy,nodes,tokensOut)

plot(phy)
phycanv <- phylocanvas(phy,treetype = "rectangular", alignlabels = F, width = 1000, textsize = 10, nodesize = 10, showcontextmenu = T, showscalebar = T)
n <- max(wt$membership)
col_vector = viridis_pal(option = "D")(n)
for (nodename in wt$names) {
  nodename1 <- gsub("\\s","_",nodename)
  #print(nodename)
  phycanv$x$nodestyles[[nodename1]]$colour <- col_vector[wt$membership[which(wt$names==nodename)]]
  #phycanv <- style_node(phycanv, nodeid = nodename, color=col_vector[wt$membership[which(wt$names==nodename)]], fillcolor=col_vector[wt$membership[which(wt$names==nodename)]])
}
phycanv

# graph to tree??
# treeLinks <- agg_links
# sources <- unique(treeLinks$Source)
# treeLinks.sub <- data.frame(Source=character(0),Target=character(0),label=character(0),weight=numeric(0),stringsAsFactors = F)
# tmp<-lapply(sources, function(x){
#   inElem <- as.data.frame(treeLinks[which(treeLinks$Source == x),])
#   inElem <- inElem[order(-inElem$weight),]
#   if(nrow(inElem) < 2){
#     treeLinks.sub <<- rbind(treeLinks.sub,inElem[1,])
#   } else{
#     treeLinks.sub <<- rbind(treeLinks.sub,inElem[1:2,])
#   }
#   
# })
# write_csv(treeLinks.sub ,"treegraph.csv")
