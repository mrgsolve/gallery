
library(minqa)
library(tidyverse)
library(mrgsolve)

data <- read_csv("indometh_data.csv", na='.')


fun <- function(p, parnames, data,pred=FALSE) {
  p <- lapply(p,exp)
  names(p) <- parnames
  mod <- param(mod, p)
  out <- mrgsim_d(mod,data, output="df", rtol=1E-12, atol=1E-12)
  if(pred) return(out)
  w <- 1/data[["DV"]]
  wres <- (data[["DV"]] - out[["CP"]])*w
  sum(wres^2,na.rm=TRUE)
  
}

mod <- modlib ("pk2")

theta <- log(c(CL = 2, V2=2, Q=2,V3=2))

fit <- newuoa(theta,fun,parnames=names(theta),data=data)

pr <- fun(fit$par,names(theta),data,pred=TRUE)

data[["PRED"]] <- pr[["CP"]]

ggplot(data, aes(PRED,DV)) + geom_point() + 
  geom_abline(intercept=0,slope=1) + 
  theme_bw() + scale_y_log10() + scale_x_log10()

library(nloptr)

system.time(
  fit <- minqa::newuoa(theta,fun,parnames=names(theta),data=data)
)
system.time(
  fitt <- nloptr::newuoa(theta,fun,parnames=names(theta),data=data)
)

system.time(
  fittt <- optim(theta,fun,parnames=names(theta),data=data)
)

system.time(
  fit <- sbplx(theta,fun,parnames=names(theta),data=data)
)

system.time(
  fit <- neldermead(theta,fun,parnames=names(theta),data=data)
)

system.time(
  fit <- auglag(theta,fun,parnames=names(theta),data=data)
)
