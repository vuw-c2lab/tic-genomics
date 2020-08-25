
runFolder <- "/Users/MLR/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/full-country-nometa_data/20/"
fils <- list.files(path = runFolder,pattern = "*.nwk")
res <- data.frame(cntry=character(0),seqs=numeric(0),diamTree=numeric(0),avgspathTree=numeric(0),diamCodon=numeric(0),avgspathCodon=numeric(0),modCodon=numeric(0),diamNucleo=numeric(0),avgspathNucleo=numeric(0),modNucleo=numeric(0))
for(fil in fils){
  tr <- read.tree(paste0(runFolder,"/",fil))
  gr <- as.igraph(tr,directed = F)
  admat <- as.matrix(read_csv(paste0(runFolder,"/",unlist(strsplit(fil,"_tree_"))[1],"_g0.csv")))
  rownames(admat) <- admat[,1]
  admat <- admat[,-1]
  gcodon <- graph_from_adjacency_matrix(admat,weighted = T)
  E(gcodon)$weight <- 1
  admat <- as.matrix(read_csv(paste0(runFolder,"/",unlist(strsplit(fil,"_tree_"))[1],"_g1.csv")))
  rownames(admat) <- admat[,1]
  admat <- admat[,-1]
  gnucleo <- graph_from_adjacency_matrix(admat,weighted = T)
  E(gnucleo)$weight <- 1
  print(unlist(strsplit(fil,"_tree_"))[1])
  print(diameter(gr))
  print(mean_distance(gr))
  res<-rbind(res,data.frame(cntry=unlist(strsplit(fil,"_tree_"))[1],seqs=length(tr$tip.label),diamTree=diameter(gr),avgspathTree=mean_distance(gr),
                            diamCodon=diameter(gcodon),avgspathCodon=mean_distance(gcodon),modCodon=modularity(gcodon,membership = cluster_walktrap(gcodon)$membership),
                            diamNucleo=diameter(gnucleo),avgspathNucleo=mean_distance(gnucleo),modNucleo=modularity(gnucleo,membership = cluster_walktrap(gnucleo)$membership)))
}
write_csv(res,"../full-country-nometa-20.csv")
tree1 <- read.tree("20_0_graphs_aligned_data/Austria_tree_nucleotides_tokens_20_0_setmetric.nwk")
tree2 <- read.tree("20_0_graphs_fulldata/Austria_tree_nucleotides_tokens_20_0_setmetric.nwk")

tree1$tip.label <- gsub("([[:digit:]]+)_-_hCoV-19/|([[:digit:]]+) - ","",tree1$tip.label)
tree2$tip.label <- gsub("([[:digit:]]+)_-_hCoV-19/|([[:digit:]]+) - ","",tree2$tip.label)
# compare trees

# using shortest paths over the common nodes
g1 <- as.igraph(tree1,directed = F)
g2 <- as.igraph(tree2, directed = F)
commonLabels <- tree2$tip.label[which(tree2$tip.label %in% tree1$tip.label)]
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


