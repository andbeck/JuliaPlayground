library(cheddar)
library(tidyverse)
library(reshape2)

#load files
sizeratioData <- read.csv("../../Results/SizeRatioTest/test.csv",header = F)

tidyr::gather(sizeratioData) %>%
  ggplot(aes(x=value,group=key,colour = key))+
    geom_line(stat = 'density')

tidyr::gather(sizeratioData) %>%
  lm(value ~ key,.) %>%
  summary()
  