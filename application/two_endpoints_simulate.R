
library(mrgsolve)
library(tidyverse)
theme_set(theme_bw())

mod <- mread_cache("two_endpoints.cpp")

set.seed(111)

out <- 
  mod %>% 
  smat(dmat(0.01,0.01)) %>%
  ev(amt = 100) %>%
  mrgsim()


plas <- 
  out %>% 
  filter(time %in% c(0.25, 0.5, 1, 2, 6, 10, 16, 24, 48)) %>%
  select(ID,time,CP) 

ur <- 
  out %>% 
  filter(time %in% c(2, 10, 20 , 30 , 40, 72)) %>% 
  select(ID,time,UR)

saveRDS(file = "two_endpoints_plasma.RDS", plas)
saveRDS(file = "two_endpoints_urine.RDS", ur)





