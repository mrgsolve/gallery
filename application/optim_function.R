
library(tidyverse)
library(mrgsolve)
theme_set(theme_bw())

source("fit_functions.R")

data <- readRDS("two_endpoints_data.RDS")

data <- filter(data, dvtype !=2)

mod <- modlib("pk1")

theta <- log(c(CL = 1, V = 10, KA = 2))

fit <- fit_newuoa(mod, dv_col="DV", pred_col="CP", theta)

ggplot() + 
  geom_point(data=fit$data, aes(time,DV)) + 
  geom_line(data=fit$data, aes(time,PRED)) + 
  scale_y_log10()



