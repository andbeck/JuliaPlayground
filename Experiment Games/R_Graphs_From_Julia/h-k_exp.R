library(tidyverse)
library(patchwork)

### h-k-t-z - APB with Eva - what is re-wiring doing.

df_r <- read_csv("outdf_r.csv")
df_nr <- read_csv("outdf_nr.csv")

glimpse(df_r)

p1 <- ggplot(df_nr, aes(x = K, y = persistence,
                colour = factor(h),
                shape = factor(h),
                group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p2 <- ggplot(df_nr, aes(x = K, y = stability,
                      colour = factor(h),
                      shape = factor(h),
                      group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p3 <- ggplot(df_nr, aes(x = K, y = diversity,
                        colour = factor(h),
                        shape = factor(h),
                        group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p4 <- ggplot(df_r, aes(x = K, y = persistence,
                        colour = factor(h),
                        shape = factor(h),
                        group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p5 <- ggplot(df_r, aes(x = K, y = stability,
                        colour = factor(h),
                        shape = factor(h),
                        group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

p6 <- ggplot(df_r, aes(x = K, y = diversity,
                        colour = factor(h),
                        shape = factor(h),
                        group = h))+
  geom_point()+
  geom_line()+
  facet_grid(`T` ~ Z)+
  scale_x_log10()

(p1+p2+p3)/(p4+p5+p6)

ggplot(df_r, aes(x = persistence, y = stability, colour = factor(h)))+
  geom_point()+
  facet_grid(Z ~ `T`)



# Masters mini project ----------
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
