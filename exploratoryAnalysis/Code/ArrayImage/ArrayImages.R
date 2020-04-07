library(cheddar)
library(tidyverse)
library(reshape2)

#load files
arrays <- list()

for(i in 1:9){
  a <- t(as.matrix(read.csv(paste0("../Results/ArrayImage/array",i,".csv"),header = F)))
  rownames(a) = paste("Species",1:20)
  colnames(a) = paste("Species",1:20)

  TL <- PredationMatrixToLinks(a)
  Sp <- data.frame(node = paste("Species",1:20))
  Wb <- list(title = i)

  arrays[[i]] <- Community(Sp,trophic.links = TL,properties = Wb)


  }

PlotWebByLevel(arrays[[1]])
PlotWebByLevel(arrays[[2]])

biomass <- read.csv("../Results/ArrayImage/Biomass.csv",header = F)
biomass$T <- 1:nrow(biomass)
biomass <- melt(biomass,id.vars = 'T')

ggplot(biomass,aes(x=T,y=value,group = variable, colour = variable))+
  geom_line()+
  theme_classic()+
  theme(legend.position = 'none')+
  xlab("Time") + ylab("Biomass")
