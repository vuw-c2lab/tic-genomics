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

meta <- read_csv("meta3.csv",col_names = T)
meta <- unique(meta)

metaTokens <- c()
for(i in 1:nrow(meta)){
  age <- ageInterval(meta$Age[i])
  district <- tolower(meta$DHB[i])
  gend <- tolower(meta$Sex[i])
  clust <- tolower(meta$outbrkno[i])
  #hosp <- meta$X14[i]
  onset <- ifelse(is.na(lubridate::mdy(meta$onsetdt[i])),toString(lubridate::as_date(meta$onsetdt[i])),toString(lubridate::mdy(meta$onsetdt[i])))
  city<- tolower(meta$lastcity[i])
  status <- tolower(meta$source[i])
  country <- tolower(meta$Travel.History[i])
  #tokenVec <- c(age,gend,clust,hosp,trans,sympt,status,country,district)
  tokenVec <- c(age,gend,clust,country,district,status,city,onset)
  #names(tokenVec) <- c("Age","Gender","Cluster","Hospitalised","Transmission","Asymptomatic","Status","LastCountry","District")
  names(tokenVec) <- c("Age","Gender","Cluster","Travel","Location","InfectionSource","LastCity","OnsetDate")
  tokenVec <- tokenVec[which(!is.na(tokenVec) & tokenVec != "NA")]
  if(length(tokenVec)>0){
    metaTokens <- c(metaTokens,paste(paste(names(tokenVec),tokenVec, sep=": "),collapse = ", "))
  } else{
    metaTokens <- c(metaTokens,"")
  }
}

meta$tokenised <- metaTokens

