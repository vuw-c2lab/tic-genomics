library(tidyverse)

ageInterval <- function(x){
  if(is.na(x)) {return(NA)}
  
  if(nchar(x)<2){
    fst <- unlist(strsplit(as.character(x), ""))[1]
    if(fst >= 5){
      return(paste0("5-9"))
    } else{
      return(paste0("0-4"))
    }
  }
  
  fst <- unlist(strsplit(as.character(x), ""))[1]
  lst <- as.numeric(unlist(strsplit(as.character(x), "")))[2]
  if(lst >= 5){
    moduA <- x %% 5
    moduB <- 10 - (x %% 10)
    return(paste0("",x-moduA,"-",x+moduB-1))
  } else{
    return(paste0("",fst,"0-",fst,"4"))
  }
}

meta <- read_csv("global_meta.csv",col_names = T)
meta$Location <- tolower(meta$Location)
meta$Gender <- tolower(meta$Gender)
meta$Specimen <- tolower(meta$Specimen)
meta$status <- tolower(meta$status)
meta$addtitional <- tolower(meta$addtitional)

metaTokens <- c()
for(i in 1:nrow(meta)){
  age <- ifelse(is.na(meta$age[i]),meta$age[i],paste0("Age: ",ageInterval(meta$age[i])))
  loc <- ifelse(is.na(meta$Location[i]),meta$Location[i],paste0("Location: ",meta$Location[i]))
  gend <- ifelse(is.na(meta$Gender[i]),meta$Gender[i],paste0("Gender: ",meta$Gender[i]))
  speci <- unlist(strsplit(meta$Specimen[i],", "))
  if(length(speci) > 1){
    speci <- paste(paste0("Specimen: ",speci),collapse=", ")
  } else if(!is.na(speci)){
    speci <- paste0("Specimen: ",speci)
  }
  #hosp <- meta$X14[i]
  #onset <- ifelse(is.na(lubridate::mdy(meta$onsetdt[i])),toString(lubridate::as_date(meta$onsetdt[i])),toString(lubridate::mdy(meta$onsetdt[i])))
  #city<- tolower(meta$lastcity[i])
  status <- unlist(strsplit(meta$status[i],", "))
  if(length(status) > 1){
    status <- paste(paste0("Status: ",status),collapse=", ")
  } else if(!is.na(status)){
    status <- paste0("Status: ",status)
  }
  addi <- unlist(strsplit(meta$addtitional[i],", "))
  if(length(addi) > 1){
    addi <- paste(paste0("Travel: ",addi),collapse=", ")
  } else if(!is.na(addi)){
    addi <- paste0("Travel: ",addi)
  }
  #tokenVec <- c(age,gend,clust,hosp,trans,sympt,status,country,district)
  tokenVec <- c(age,gend,loc,speci,status,addi)
  #names(tokenVec) <- c("Age","Gender","Cluster","Hospitalised","Transmission","Asymptomatic","Status","LastCountry","District")
  names(tokenVec) <- c("Age","Gender","Location","Specimen","Status","Travel")
  tokenVec <- tokenVec[which(!is.na(tokenVec) & tokenVec != "NA")]
  if(length(tokenVec)>0){
    metaTokens <- c(metaTokens,paste(tokenVec,collapse = ", "))
  } else{
    metaTokens <- c(metaTokens,"")
  }
}

meta$tokenised <- metaTokens





#### NZ META


metaNZ <- read_csv("meta3.csv",col_names = T)
metaNZ <- unique(metaNZ)

metaTokens <- c()
for(i in 1:nrow(metaNZ)){
  age <- ageInterval(metaNZ$Age[i])
  loc <- metaNZ$DHB[i]
  gend <- metaNZ$Sex[i]
  clust <- metaNZ$outbrkno[i]
  #hosp <- meta$X14[i]
  #onset <- ifelse(is.na(lubridate::mdy(meta$onsetdt[i])),toString(lubridate::as_date(meta$onsetdt[i])),toString(lubridate::mdy(meta$onsetdt[i])))
  city<- tolower(metaNZ$lastcity[i])
  status <- metaNZ$source[i]
  #country <- tolower(metaNZ$Travel.History[i])
  #tokenVec <- c(age,gend,clust,hosp,trans,sympt,status,country,district)
  tokenVec <- c(age,gend,city,loc)
  #names(tokenVec) <- c("Age","Gender","Cluster","Hospitalised","Transmission","Asymptomatic","Status","LastCountry","District")
  names(tokenVec) <- c("Age","Gender","Travel","Location")
  tokenVec <- tokenVec[which(!is.na(tokenVec) & tokenVec != "NA")]
  if(length(tokenVec)>0){
    metaTokens <- c(metaTokens,paste(paste(names(tokenVec),tokenVec, sep=": "),collapse = ", "))
  } else{
    metaTokens <- c(metaTokens,"")
  }
}

metaNZ$tokenised <- metaTokens

