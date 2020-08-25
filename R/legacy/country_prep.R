library("Biostrings")
library(stringr)
library(lubridate)
library(tesseract)
library(readr)
library(msa)
library(DECIPHER)
setwd("~/OneDrive - Victoria University of Wellington - STAFF/RScripts/2020-03-09-epidemic-dashboard/")

ds <- readRDS("global_nometa_aligned.rda")
subset.ds <- readRDS("global_nometa_aligned_subset.rda")

rankedCountries <- plyr::count(unlist(lapply(strsplit(ds$nams,"/"),function(x){x[2]})))
rankedCountries <- rankedCountries[order(-rankedCountries$freq),]

cutoff <- min(rankedCountries$freq[c(10,14,15,16,23,26,29,35)])

rankedCountries <- as.character(rankedCountries$x[c(10,14,15,16,23,26,29,35)])
#rankedCountries <- c(rankedCountries,"New Zealand")
#which(!grepl(paste(rankedCountries,collapse="|"),subset.ds$nams))
#subset.ds <- subset.ds[which(!grepl(paste(rankedCountries,collapse="|"),subset.ds$nams)),]

for(country in rankedCountries){
  ds <- readRDS("global_nometa_aligned.rda")
  ds <- ds[which(grepl(paste0(country),ds$nams)),]
  #ds <- ds[1:cutoff,]
  
  #for(i in 1:nrow(subset.ds)){
  #  ds <- rbind(ds,data.frame(seq=subset.ds$seq[i],dats=subset.ds$dats[i],nams=subset.ds$nams[i],metatoks=subset.ds$metatoks[i]))
  #}
  ds$dats <- as.character(ds$dats)
  ds$nams <- as.character(ds$nams)
  ds <- ds[order(ds$dats),]
  print(paste0(country,nrow(ds)))
  ds_un<-unique(data.table::setDT(ds), by = c("nams"))
  saveRDS(ds,paste0(gsub(" ","_",tools::toTitleCase(country)),"_aligned.rda"))
  names <- ds$nams
  
  write_csv(as.data.frame(ds$seq),paste0("ordered_",gsub(" ","_",tools::toTitleCase(country)),".txt"), col_names = F)
  write_csv(as.data.frame(ds$metatoks),paste0("ordered_",gsub(" ","_",tools::toTitleCase(country)),"_meta.txt"), col_names = F)
  saveRDS(as.data.frame(ds$seq),paste0(gsub(" ","_",tools::toTitleCase(country)),"-seq.rda"))
}
