library(ggplot2)
library(dplyr)

d <- data.frame (
  Training = c("Strength", "Stamina", "Other"),
  Pulse = c(100, 150, 120),
  Duration = c(60, 30, 45)
)

d %>% ggplot( aes(x=Pulse, y=Duration)) +
  geom_line() +
  geom_point()