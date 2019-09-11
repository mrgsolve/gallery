Sobol sensitivity analysis
================
Kyle Baron and Ahmed Elmokadem
2019-09-04 14:10:20

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

    .          TVCL       TVVC         TVKA         TVQ       TVVP
    . 1    1.361773  559.02011  1.411879296   6.5166076   12.45871
    . 2    7.528564 5791.84475  0.003626071 137.5590446 8153.91343
    . 3    2.838735 6238.19311  0.006786784   3.3999461   40.93909
    . 4 1083.214399  596.20614 10.351813055  27.3969644  257.47784
    . 5  222.154810 7873.06062  0.025234897   0.5845153  100.41110
    . 6   29.087651   27.43569  0.004670174   0.1685278 2471.69992

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
    .         original          bias  std. error   min. c.i.  max. c.i.
    . TVCL 0.111610862 -1.718937e-04 0.006474346 0.098965610 0.12602305
    . TVVC 0.153619118 -9.389818e-04 0.008221707 0.139729172 0.17155083
    . TVKA 0.028567213  1.872252e-04 0.002373118 0.023773938 0.03313045
    . TVQ  0.008698163  2.222233e-05 0.002046892 0.004592815 0.01244719
    . TVVP 0.004631098 -6.015749e-05 0.001256271 0.002085233 0.00704924
    . 
    . Total indices:
    .        original          bias  std. error  min. c.i. max. c.i.
    . TVCL 0.69602357 -3.135509e-04 0.011714859 0.67593136 0.7220769
    . TVVC 0.74711857 -8.026551e-04 0.010320519 0.72897110 0.7679636
    . TVKA 0.25952423 -2.597943e-05 0.011324677 0.23568274 0.2840098
    . TVQ  0.22073262  2.781637e-04 0.009763087 0.20292239 0.2407872
    . TVVP 0.09843795  1.719443e-04 0.007096819 0.08291732 0.1120978
