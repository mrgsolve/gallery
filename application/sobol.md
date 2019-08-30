Sobol sensitivity analysis
================
Kyle Baron and Ahmed Elmokadem
2019-08-30 07:42:52

  - [Reference / About](#reference-about)
  - [Tools](#tools)
  - [The sunitinib PK model](#the-sunitinib-pk-model)
      - [Sunitinib dosing](#sunitinib-dosing)
  - [Generate samples](#generate-samples)
  - [A bunch of helper functions](#a-bunch-of-helper-functions)
  - [Run the analysis](#run-the-analysis)
      - [First, generate the samples](#first-generate-the-samples)
      - [Then, run
        `sensitivity::sobol2007`](#then-run-sensitivitysobol2007)
  - [Results](#results)

# Reference / About

Zhang XY, Trame MN, Lesko LJ, Schmidt S. **Sobol Sensitivity Analysis: A
Tool to Guide the Development and Evaluation of Systems Pharmacology
Models**. CPT Pharmacometrics Syst Pharmacol. 2015 Feb;4(2):69-79. doi:
10.1002/psp4.6. PubMed PMID:
[27548289](https://www.ncbi.nlm.nih.gov/pubmed/27548289)

This example replicates an analysis presented in the Zhang et al. paper,
but here using mrgsolve and other tools available for R.

# Tools

``` r
library(mrgsolve)
library(tidyverse)
library(PKPDmisc)
library(sensitivity)
```

# The sunitinib PK model

``` r
mod <- mread("sunit", "model") %>% 
  update(end = 24, delta = 1) %>% zero_re
```

``` r
see(mod)
```

    . 
    . Model file:  sunit.cpp 
    . $PARAM
    . TVCL = 51.8
    . TVVC = 2030
    . TVKA = 0.195
    . TVQ = 7.22
    . TVVP = 583
    . WTVC = 0.459
    . SEXCL = -0.0876
    . ASIANCL = -0.130
    . GISTCL = -0.285
    . SOLIDCL = -0.269
    . MRCCCL = -0.258
    . SEX = 0, ASIAN = 0, GIST = 0
    . SOLID = 0, MRCC = 0, WT = 76.9
    . 
    . $MAIN
    . double CL  = TVCL * (1+SEXCL*SEX) * (1+ASIANCL*ASIAN) * 
    .   (1+GISTCL*GIST) * (1+SOLIDCL*SOLID) * (1+MRCCCL*MRCC) * exp(ETA(1));
    . 
    . double V2 = TVVC*pow(WT/76.9, WTVC)*exp(ETA(2));
    . double KA = TVKA*exp(ETA(3));
    . double Q  = TVQ;
    . double V3 = TVVP;
    . 
    . $OMEGA 0.14 0.18 0.64
    . 
    . $SIGMA 0.146
    . 
    . $PKMODEL cmt = "GUT CENT, PERIPH", depot = TRUE
    . 
    . $POST
    . capture CP = (1000*CENT/V2);

## Sunitinib dosing

``` r
sunev <- function(amt = 50,...) ev(amt = amt, ...)
```

# Generate samples

Th function generates uniform samples from a 100 fold decrease to 100
fold increase in the nominal parameter value.

The return value is a list with two data frames that can be passed into
the sobol function.

``` r
gen_samples <- function(n, l, which = names(l), 
                        factor = c(0.01,100)) {
  
  vars <- tidyselect::vars_select(names(l), !!(enquo(which)))
  
  l <- as.list(l)[vars]
  
  l <- lapply(l, function(x) {
    x*factor  
  })
  
  n <- length(l)*n*2
  
  df <- as.data.frame(l)
  
  len <- length(df)
  
  X <- matrix(ncol=len, nrow=n)
  
  colnames(X) <- names(df)
  
  Y <- X
  
  for(i in seq(len)){
    r <- runif(n, df[1,i], df[2,i])
    X[,i] <- r
    r <- runif(n, df[1,i], df[2,i])
    Y[,i] <- r
  }
  
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}
```

# A bunch of helper functions

Simulate one chunk of data; not absolutely needed for this; but we might
use it if we expand this to run in parallel. Ill keep it around for now.

``` r
sim_chunk <- function(mod, x) {
  mrgsim_ei(x = mod, ev = sunev(), idata = x, obsonly = TRUE, output="df")
}
```

Simulate a batch of data. The summary is AUC for each individual.

``` r
batch_run <- function(x) {
  out <- sim_chunk(mod,x)
  out <- 
    group_by(out,ID) %>% 
    summarise(AUC = auc_partial(time,CP))
  return(out$AUC)
}
```

# Run the analysis

## First, generate the samples

``` r
samp <- gen_samples(10000, param(mod), TVCL:TVVP)
head(samp$x1)
```

    .        TVCL       TVVC      TVKA       TVQ      TVVP
    . 1  995.1216   8088.664 15.189657 528.73175  9207.374
    . 2 2143.7808 104116.499  7.707399  13.02232  9378.689
    . 3 2878.4971  14088.020  2.606306 339.97947 53009.154
    . 4  859.8960  86158.704  7.157743 234.57548 44027.799
    . 5 4435.2362  30600.809 17.004325 464.45263 40380.444
    . 6 4971.4705  10292.308 14.653676 224.10473 24536.750

## Then, run `sensitivity::sobol2007`

``` r
x <- sobol2007(batch_run, X1=samp$x1, X2=samp$x2, nboot=100)
```

# Results

``` r
plot(x)
```

![](img/sobolunnamed-chunk-11-1.png)<!-- -->

``` r
x
```

    . 
    . Call:
    . sobol2007(model = batch_run, X1 = samp$x1, X2 = samp$x2, nboot = 100)
    . 
    . Model runs: 700000 
    . 
    . First order indices:
    .           original          bias   std. error    min. c.i.   max. c.i.
    . TVCL  0.1438667922  5.341818e-04 0.0119131431  0.120037569 0.166384912
    . TVVC  0.3313325821  7.191619e-04 0.0286001529  0.273259327 0.385705047
    . TVKA  0.0022658887 -9.265332e-05 0.0004553831  0.001405788 0.003191302
    . TVQ   0.0015259725  6.506356e-05 0.0015534711 -0.001800230 0.004651025
    . TVVP -0.0008613722  2.024928e-04 0.0009967453 -0.002930682 0.001200853
    . 
    . Total indices:
    .         original          bias  std. error    min. c.i.  max. c.i.
    . TVCL 0.673032352  0.0013062321 0.022961641 6.226432e-01 0.71598341
    . TVVC 0.838137911  0.0004124743 0.012924918 8.117653e-01 0.86460619
    . TVKA 0.006402944  0.0001057573 0.002785053 9.205852e-05 0.01141503
    . TVQ  0.096955907  0.0003735636 0.041281016 5.598757e-03 0.16590543
    . TVVP 0.032099945 -0.0005689257 0.014729286 3.378184e-03 0.05908883
