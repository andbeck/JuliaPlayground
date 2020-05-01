library(ggplot2)
library(readr)
library(reshape2)

test <- read.csv("../../Results/FirstTest/testGilljam.csv",header = F)

colnames(test) <- c("NoRewiring","GiljamRewiring","ADBMRewiring")

a <- melt(test)

plt <- ggplot(a,aes(x = value,group = variable,colour = variable))+
  geom_line(stat = 'density')+
  xlab("Proportion of Remaining Species")+
  theme_classic()

ggsave("../../Results/FirstTest/test.pdf",plot = plt)
write.csv(test,"../../Results/FirstTest/test.csv")
