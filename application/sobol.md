Sobol sensitivity analysis
================
Kyle Baron and Ahmed Elmokadem
2021-07-28 17:35:08

-   [Reference / About](#reference--about)
-   [Tools](#tools)
-   [The sunitinib PK model](#the-sunitinib-pk-model)
    -   [Sunitinib dosing](#sunitinib-dosing)
-   [Generate samples](#generate-samples)
-   [A bunch of helper functions](#a-bunch-of-helper-functions)
-   [Run the analysis](#run-the-analysis)
    -   [First, generate the samples](#first-generate-the-samples)
    -   [Then, run
        `sensitivity::sobol2007`](#then-run-sensitivitysobol2007)
-   [Results](#results)

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

Th function generates uniform samples from a 5 fold decrease to 5 fold
increase in the nominal parameter value.

The return value is a list with two data frames that can be passed into
the sobol function.

``` r
gen_samples <- function(n, l, which = names(l), 
                        factor = c(0.1,10)) {
  
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
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    X[,i] <- r
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
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

    .         TVCL       TVVC       TVKA       TVQ       TVVP
    . 1   8.608178  3319.7761 0.02294751  2.441683  120.16921
    . 2   8.713837 15791.3948 0.63534638 10.709428 1168.69128
    . 3 368.959778  1056.8573 0.02819864 24.773588   97.78583
    . 4 400.369944  1376.2384 0.13515809  7.143678  818.15172
    . 5 341.677497   459.6150 0.18384088 19.358219  562.02095
    . 6 181.220458   673.0695 0.03846543  7.136121 2476.55036

## Then, run `sensitivity::sobol2007`

``` r
x <- sobol2007(batch_run, X1=samp$x1, X2=samp$x2, nboot=1000)
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
    . sobol2007(model = batch_run, X1 = samp$x1, X2 = samp$x2, nboot = 1000)
    . 
    . Model runs: 700000 
    . 
    . First order indices:
    .          original          bias   std. error   min. c.i.   max. c.i.
    . TVCL 0.1979690453 -1.053888e-04 0.0052982790  0.18764614 0.208241256
    . TVVC 0.3750716575  2.058715e-05 0.0080033687  0.35936809 0.391380259
    . TVKA 0.0669549268  8.332064e-05 0.0023052791  0.06248294 0.071461807
    . TVQ  0.0057798916  1.706376e-05 0.0010712851  0.00366490 0.007815363
    . TVVP 0.0002899744  6.098216e-06 0.0006135896 -0.00098076 0.001465525
    . 
    . Total indices:
    .        original          bias  std. error   min. c.i.  max. c.i.
    . TVCL 0.49950899 -3.385304e-05 0.007106738 0.486126040 0.51325048
    . TVVC 0.70335670 -2.747905e-05 0.006708954 0.689897688 0.71692880
    . TVKA 0.15247316 -4.963958e-05 0.005027915 0.142499870 0.16235736
    . TVQ  0.03812596 -1.121061e-04 0.003581700 0.031257933 0.04537120
    . TVVP 0.01225372 -8.378933e-05 0.001861244 0.008611578 0.01588701
