library(digest)
library(entropy)
library(vegan)
createFeatures <- function(nodes,type=NULL){
  #shannon_entropy(Dist(plyr::count(foo[[1]])$freq))
  #allTokens<-unique(unlist(str_split_fixed(nodes$token,", ",n = Inf)[1,]))
  #inter <- setNames(as.list(rep(0,length(allTokens))),allTokens)
  #wien <- data.frame(ShannonWiener=numeric(nrow(nodes)),Pielou=numeric(nrow(nodes)),Richness=numeric(nrow(nodes)))
  #colnames(links) <- c('source', 'target', 'tag', 'weight')
  allTokens<-unique(unlist(strsplit(nodes$token,", ")))
  if(type=="metadata"){
    allTokens <- allTokens[which(grepl("Age: ",allTokens) | grepl("Gender: ",allTokens) | grepl("Location: ",allTokens) | grepl("Travel: ",allTokens))]
  }
  inter <- setNames(as.list(rep(0,length(allTokens))),allTokens)
  coordinates <- data.frame(t=numeric(nrow(nodes)),specificity=numeric(nrow(nodes)),diversity=numeric(nrow(nodes)))
  wien <- data.frame(ShannonWiener=numeric(nrow(nodes)),Pielou=numeric(nrow(nodes)),Richness=numeric(nrow(nodes)))
  spec <- list()
  div <- 1
  
  #z <- 0
  
  res <- unlist(lapply(1:nrow(nodes), function(x,y){
    print("####")
    t1<-Sys.time()
    print(t1)
    if(type=="metadata"){
      checkMeta <- grepl("Age: ",nodes[x,2]) | grepl("Gender: ",nodes[x,2]) | grepl("Location: ",nodes[x,2]) | grepl("Travel: ",nodes[x,2])
    } else{
      checkMeta <- T
    }
    if(nodes[x,2] == "" | is.na(nodes[x,2]) | !checkMeta){
      coordinates[x,] <<- c(x, 0, 0)
      if(x > 1){
        #ent <- rbind(ent, ent[nrow(ent),])
        wien[x,] <<- wien[x-1,]

      }else{
        #ent <- rbind(ent, c(0, 1, 0, 1))
        wien[x,] <<- c(0, 1, 1)
      }
    } else{
      theInter <- sort(unlist(strsplit(as.character(nodes[x,2]), ', ')))
      if(type=="metadata"){
        theInter <- theInter[which(grepl("Age: ",allTokens) | grepl("Gender: ",allTokens) | grepl("Location: ",allTokens) | grepl("Travel: ",allTokens))]
      }
      sortedInter <- paste(theInter, collapse = ', ')
      nextI <- digest::digest(sortedInter, algo = "md5")
      if(length(spec) == 0){
        coordinates[x,] <<- c(x, 1, 1)
        spec[[nextI]] <<- c(1, 1)
      } else{
        if(is.null(spec[[nextI]])){
          div <<- div+1
          spec[[nextI]] <<- c(1,div)
          coordinates[x,] <<- c(x,spec[[nextI]][1],div)
        }else{
          spec[[nextI]] <<- c(spec[[nextI]][1]+1,spec[[nextI]][2])
          coordinates[x,] <<- c(x,spec[[nextI]][1],spec[[nextI]][2])
        }
      }
      print(Sys.time()-t1)
      interact <- list()
      cooccure <- list()
      #temp1 <- c()
      interactions <- unlist(strsplit(unlist(sortedInter),', '))
      for(v in 1:length(interactions)){
        inter[[interactions[v]]] <<- inter[[interactions[v]]] + 1
      }
      print(length(inter[inter==0]))
      df <- data.frame(unlist(inter[inter!=0]))
      print(nrow(df))
      # tmp <- df[,1] / colSums(df)
      # df$loga <- log(tmp)
      # df$piloga <- tmp * log(tmp)
      # if(is.nan((-1 * (colSums(df)[3])) / log(nrow(df)))){
      #   ent <- rbind(ent,c(entropy.empirical(df[,1], unit = "log2"), 1, -1 * (colSums(df)[3]), 1))
      # } else{
      #   ent <- rbind(ent, c(entropy.empirical(df[,1], unit = "log2"), (-1 * (colSums(df)[3])) / log2(nrow(df)),-1*(colSums(df)[3]),(-1*(colSums(df)[3]))/log(nrow(df))))
      # }

      if(nrow(df) == 1){
        wien[x,] <<- c(0, 1, 1)
      } else{
        S <- nrow(df)
        H <- vegan::diversity(df[,1])
        J <- H/log(S)
        wien[x,] <<- c(H, J, S)
      }
      
    }
    #distr <- rinform::shannon_entropy(rinform::Dist(plyr::count(unlist(str_split(nodes[1:x,2],", ", simplify = F)))$freq))
    print(Sys.time()-t1)
    #distr
  },y=nodes))
  
  coordinates <- cbind(coordinates,wien)
  
  return(coordinates)
}

setwd("/Users/MLR/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/full-country-nometa_data/")

folders <- c("New_Zealand")
folders <- c("Saudi_Arabia","Austria","France","Belgium","Japan","Colombia","Russia","Taiwan","New_Zealand")
tokenisers <- c("codon","nucleotides")
for(i in folders){
  for(j in tokenisers){
    #1st layer
    links <- read_csv(paste0(i,"/",j,"/links.csv"))
    nodes <- read_csv(paste0(i,"/",j,"/nodes.csv"))
    
    colnames(links) <- c("id1","id2","label")
    
    ds <- readRDS(paste0(i,"_aligned.rda"))
    names <- ds$nams
    colnames(nodes) <- c("id","token","ord")
    nodes$genenames <- paste(nodes$id," - ",names,sep="")
    
    #saveRDS(links,paste0(tools::toTitleCase(i),"-",j,"-links.rda"))
    #saveRDS(nodes,paste0(tools::toTitleCase(i),"-",j,"-nodes.rda"))
    if(j=="metadata"){
      features <- createFeatures(nodes,type = "metadata")
    } else{
      features <- createFeatures(nodes, type = "genomes")
    }
    saveRDS(features,paste0(tools::toTitleCase(i),"-",j,"-features.rda"))
  }
}
