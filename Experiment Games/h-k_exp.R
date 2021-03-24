library(tidyverse)
library(patchwork)

df <- read_csv("h-k_exp.csv")

sumDat_div <- df %>% 
  group_by(h,K) %>% 
  summarise(
    meanDiv = mean(diversity, na.rm = TRUE),
    meanBiomass = mean(biomass, na.rm = TRUE),
    meanStability = mean(stability, na.rm = TRUE),
  )
  
p1 <- ggplot(sumDat_div, aes(x = K, y = meanDiv, colour = factor(h)))+
  geom_line()+
  scale_x_log10()

p2 <- ggplot(sumDat_div, aes(x = K, y = meanBiomass, colour = factor(h)))+
  geom_line()+
  scale_x_log10()

p3 <- ggplot(sumDat_div, aes(x = K, y = meanStability, colour = factor(h)))+
  geom_line()+
  scale_x_log10()

p1+p2+p3

### h-k-t-z

df2 <- read_csv("outdf.csv")
glimpse(df2)
df2 %>% distinct(`T`)
p1 <- ggplot(df2, aes(x = K, y = persistence,
                colour = factor(h),
                shape = factor(h),
                group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p2 <- ggplot(df2, aes(x = K, y = stability,
                      colour = factor(h),
                      shape = factor(h),
                      group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p1+p2

ggplot(df2, aes(x = persistence, y = stability, colour = factor(h)))+
  geom_point()+
  facet_grid(K ~ `T`+Z)
