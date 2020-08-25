library("Biostrings")
library(stringr)
library(lubridate)
library(tesseract)
library(readr)
library(DECIPHER)
library(e1071)
set.seed(3011)
setwd("~/OneDrive - Victoria University of Wellington - STAFF/RScripts/2020-03-09-epidemic-dashboard/")
ref_genome <- as.character(read_csv("ref_genome.txt",col_names = F)[,1])

#set <- readDNAStringSet("gisaid_hcov-19_2020_05_20_12_global_with_meta.fasta")
# high coverage only
#set <- readDNAStringSet("gisaid_hcov-19_2020_05_11_01.fasta")
# low coverage excluded
sets <- readDNAStringSet("NZCovid/nz_data/global/aligned.fasta")
meta <- read_csv("NZCovid/nz_data/global/metadata.csv",col_names = T)

whereFrom <- c("Brazil","England","USA","Australia","Belgium","Spain","China","Italy","NewZealand","Netherlands","Poland","India","HongKong","UnitedArabEmirates","Norway","Taiwan","Pakistan","Austria","Sweden","Scotland","Wales","France")

seqs <- names(sets)


#ds <- data.frame(seq=seqs[which(gsub("hCoV-19/","",names) %in% nxtree$tip.label[which(grepl("/",nxtree$tip.label))])],dats=dats[which(gsub("hCoV-19/","",names) %in% nxtree$tip.label[which(grepl("/",nxtree$tip.label))])],nams=names[which(gsub("hCoV-19/","",names) %in% nxtree$tip.label[which(grepl("/",nxtree$tip.label))])],stringsAsFactors = F)
prep <- data.frame(seq=character(0),dats=character(0),nams=character(0),stringsAsFactors = F)

yVec <- unlist(strsplit(substr(ref_genome,55,29797),""))

for(idx in 1:length(sets)){
  set <- sets[idx]
  nam <- seqs[idx]
  seq <- substr(toString(set),55,29797)
  
  if(nchar(seq) == length(yVec)){
  
    xVec <- unlist(strsplit(seq,""))
    dist <- hamming.distance(xVec,yVec)
    
    match <- regexpr("(N)\\1+",seq)
    print(attributes(match)$match.length)
    #if(dist<=5000 & max(attributes(match)$match.length < 100)){
      getMeta <- meta[which(meta$strain==nam),]
      if(grepl("/",getMeta$date) & (unlist(strsplit(nam,"/"))[1] %in% whereFrom)){
        #if(lubridate::dmy(getMeta$date) > "2020-04-01"){
          prep <- rbind(prep,data.frame(seq=seq,dats=lubridate::dmy(getMeta$date),nams=nam))
        #}
      }
    #}
  }
}
prep$dats <- as.character(prep$dats)
prep$nams <- as.character(prep$nams)
prep <- prep[order(prep$dats),]

saveRDS(prep,"plain_global_aligned.rda")

