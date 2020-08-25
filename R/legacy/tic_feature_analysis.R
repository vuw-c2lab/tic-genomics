library(tidyverse)
library(plotly)
library(ggthemes)
library(ggpmisc)
library(cowplot)
setwd("/Users/MLR/OneDrive - Victoria University of Wellington - STAFF/Python Scripts/tic/output/aligned-country_data")
folders <- list("Saudi_Arabia","Austria","France","Belgium","Japan","Colombia","Russia","Taiwan","New_Zealand")
tokenisers <- list("codon","nucleotides","metadata")

fetOut <- NULL

for(i in folders){
  for(j in tokenisers){
    fet <- readRDS(paste0(i,"-",j,"-features.rda"))
    alin <- readRDS(paste0(i,"_aligned.rda"))
    fet$country <- i
    fet$tokeniser <- j
    if(is.null(alin$ord)){
      fet$dtime <- alin$dats
    } else{
      fet$dtime <- alin$ord
    }
    
    if(is.null(fetOut)) {
      fetOut <- fet
    } else{
      fetOut <- rbind(fetOut,fet)
    }
    # plot_ly(x=fet$t,
    #         y=fet$ShannonWiener,
    #         z = fet$Richness, 
    #         type = "scatter3d",mode="points",opacity=0.5)
  }
}
peaksandvals <- NULL
plotList <- list()
for(i in folders){
  a <- length(quantmod::findValleys(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country==i),4],thresh = .00001))
  b <- length(quantmod::findPeaks(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country==i),4],thresh = .00001))
  if(is.null(peaksandvals)){
    peaksandvals <- data.frame(country=i,peaks=b,valleys=a)
  } else{
    peaksandvals <- rbind(peaksandvals,data.frame(country=i,peaks=b,valleys=a))
  }
  
  plotDtime <- fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country==i),]
  plotDtime$dtime <- as.Date(plotDtime$dtime)
  data.length <- length(plotDtime$dtime)
  time.min <- plotDtime$dtime[1]
  time.max <- plotDtime$dtime[data.length]
  all.dates <- seq(time.min, time.max, by="day")
  all.dates.frame <- data.frame(list(time=all.dates))
  merged.data <- merge(all.dates.frame, plotDtime, all=T)
  merged.data$ShannonWiener[which(is.na(merged.data$ShannonWiener))] <- 0
  
  plotList[[i]] <- ggplot(merged.data,aes(x=dtime,y=ShannonWiener)) +
    geom_point() + theme_clean() +
    stat_smooth(method = lm, formula = y ~ splines::bs(x, df = 3)) +
    ylab("Entropy") + xlab("Time")
}

plot_grid(plotList[[1]], plotList[[2]], 
          plotList[[3]], plotList[[4]],
          plotList[[5]], plotList[[6]],
          plotList[[7]], plotList[[8]],plotList[[9]],ncol = 3,labels=names(plotList),
          label_size = 10,
          label_x = 0, label_y = 0,
          hjust = -0.5, vjust = -0.5)

wilcox.test(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="New_Zealand"),4],fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Austria"),4])
ks.test(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="New_Zealand"),4],fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Austria"),4])
#kruskal.test(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="New_Zealand"),4],fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Taiwan"),4])

#### BACKUP ----
cutoff <- data.frame( x = c(-Inf, Inf), y = max(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Belgium"),4]), cutoff = factor(max(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Belgium"),4])) )
ggplot(fetOut[which(fetOut$tokeniser=="nucleotides"),],aes(x=t,y=Pielou,colour=country)) +
  geom_point() + theme_clean() +
  stat_smooth(method = lm, formula = y ~ splines::bs(x, df = 5)) +
  #geom_line(aes(x, y), cutoff, linetype="dashed") +
  theme(legend.position="none")


pop.ss <- nls(ShannonWiener ~ SSlogis(c(1:nrow(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Saudi_Arabia"),])), phi1, phi2, phi3), data = fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Saudi_Arabia"),])
predict.pop.ss <- predict(pop.ss, c(1:nrow(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Saudi_Arabia"),])))
ggplot(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Saudi_Arabia"),],
       aes(x=t,y=ShannonWiener,colour=country)) +
  geom_line() +
  theme_clean() + stat_peaks(colour = "red",span = 9) + geom_line(aes(y=predict.pop.ss))

pracma::findpeaks(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="New_Zealand"),4])

quantmod::findPeaks(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="New_Zealand"),4],thresh = .000001)

quantmod::findValleys(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="Taiwan"),4],thresh = .000001)

splus2R::peaks(fetOut[which(fetOut$tokeniser=="nucleotides" & fetOut$country=="New_Zealand"),4],span = 5)