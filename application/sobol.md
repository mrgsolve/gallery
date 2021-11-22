Sobol sensitivity analysis
================
Kyle Baron and Ahmed Elmokadem
2021-11-22 16:01:02

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

**NOTE**: This example uses the `sensitivity` package to run the
analysis. This package works well, but I am now preferring the
`sensobol` package. You can see the same analysis run with `sensobol`
here: [global-sensobol.md](global-sensobol.md).

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
                        factor = c(0.2,5)) {
  
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

    .        TVCL      TVVC       TVKA      TVQ      TVVP
    . 1 105.19835 3683.1006 0.04088887 2.334924  549.2293
    . 2  82.81668 3745.2574 0.04704564 2.985555  156.0289
    . 3  80.44561  564.3977 0.30977500 5.207510  517.3773
    . 4  40.43505 8619.9091 0.19579378 2.128595  169.9813
    . 5  76.97500 3152.2796 0.34911939 1.901407 1803.5818
    . 6  99.79838  629.4909 0.07681993 1.996189  397.5637

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
    .          original         bias   std. error    min. c.i.   max. c.i.
    . TVCL 0.1823489631 2.657960e-05 0.0050160337 1.725997e-01 0.192263576
    . TVVC 0.5141753555 3.220217e-04 0.0084467828 4.975119e-01 0.532512447
    . TVKA 0.0749283884 2.617883e-05 0.0023733276 7.011287e-02 0.079424818
    . TVQ  0.0035332569 3.077143e-05 0.0007774327 2.012456e-03 0.005085881
    . TVVP 0.0007372073 7.951747e-06 0.0003337387 6.274435e-05 0.001391284
    . 
    . Total indices:
    .         original          bias   std. error    min. c.i.   max. c.i.
    . TVCL 0.377112176 -8.662422e-05 0.0063352451  0.365162809 0.389724738
    . TVVC 0.720157338 -5.407326e-05 0.0069729474  0.705741837 0.733569927
    . TVKA 0.112670286  8.308806e-05 0.0040604439  0.104782560 0.120670548
    . TVQ  0.007889849 -6.996931e-05 0.0018116421  0.004318470 0.011529577
    . TVVP 0.001242264 -5.662055e-06 0.0007009531 -0.000142662 0.002668531
