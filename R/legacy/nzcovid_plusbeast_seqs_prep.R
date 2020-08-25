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

subset.ds <- readRDS("global_aligned_beast_seqs.rda")

fils <- list.files("NZCovid/nz_data/genomes",full.names = T)

prep <- data.frame(seq=character(0),ord=numeric(0),nam=character(0),stringsAsFactors = F)

names <- data.frame(nam=character(0),stringsAsFactors = F)

#meta <- read_csv("meta.csv",col_names = F)
yVec <- unlist(strsplit(substr(ref_genome,55,29797),""))

for(fil in fils){
  set <- readDNAStringSet(fil)
  nam <- gsub("NZCovid/nz_data/genomes/","",gsub(".consensus.fasta","",fil))
  seq <- substr(toString(set),55,29797)
  xVec <- unlist(strsplit(seq,""))
  dist <- hamming.distance(xVec,yVec)
  
  if(nam %in% nxtree$tip.label[which(!grepl("/",nxtree$tip.label))]){
    if(length(which(meta$Accession==nam))>0){
      if(!is.na(meta$onsetdt[which(meta$Accession==nam)]) | !is.na(meta$Date.Collected[which(meta$Accession==nam)])){
        if(!is.na(meta$onsetdt[which(meta$Accession==nam)])){
          dat <- ifelse(is.na(lubridate::mdy(meta$onsetdt[which(meta$Accession==nam)])),toString(lubridate::ymd(meta$onsetdt[which(meta$Accession==nam)])),toString(lubridate::mdy(meta$onsetdt[which(meta$Accession==nam)])))
        } else {
          dat <- toString(meta$Date.Collected[which(meta$Accession==nam)])
        }
        prep <- rbind(prep,data.frame(seq=seq,ord=dat,nam=nam))
      }
    }
  }
  #}
}

prep <- prep[which(!is.na(prep$ord)),]


for(i in 1:nrow(subset.ds)){
  #if(!grepl("New Zealand",subset.ds$nams[i]) & !is.na(subset.ds$nams[i])){
    prep <- rbind(prep,data.frame(seq=subset.ds$seq[i],ord=subset.ds$dats[i],nam=subset.ds$nams[i]))
  #}
}

prep$nam <- as.character(prep$nam)
prep$seq <- as.character(prep$seq)
prep$ord <- as.character(prep$ord)
#prep$metatoks <- as.character(prep$metatoks)

prep <- unique(data.table::setDT(prep), by = c("nam"))

prep <- prep[order(prep$ord),]
names <- prep$nam


names <- as.character(names)

saveRDS(as.data.frame(prep$seq),"New_Zealand2-seq.rda")
saveRDS(prep,"New_Zealand2_aligned.rda")

write_csv(as.data.frame(prep$seq),"ordered_nz_plus_beast.txt",col_names = F)
