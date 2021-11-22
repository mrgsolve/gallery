Sobol sensitivity analysis using sensobol
================
Kyle Baron
2021-11-22 15:49:08

-   [Tools](#tools)
-   [The sunitinib PK model](#the-sunitinib-pk-model)
    -   [Sunitinib dosing](#sunitinib-dosing)
-   [Generate samples](#generate-samples)
-   [Run the analysis](#run-the-analysis)
    -   [Simulation](#simulation)
    -   [Calculate AUC](#calculate-auc)
    -   [Indices](#indices)
    -   [Visualize](#visualize)

# Tools

``` r
library(mrgsolve)
library(tidyverse)
library(PKPDmisc)
library(data.table)
library(sensobol)
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

**Number**

``` r
N <- 2 ^ 15
```

**Distribution**

For this example, we will just take uniform samples based on the
logarithm of the current value of the parameter.

``` r
params <- unlist(as.list(param(mod))[c("TVCL", "TVVC", "TVKA", "TVQ", "TVVP")])
umin <- log(params / 5)
umax <- log(params * 5)
```

**Generate**

``` r
mat <- as_tibble(sobol_matrices(N = N, params = c("TVCL", "TVVC", "TVKA", "TVQ", "TVVP")))
```

**Transform and groom**

``` r
mat2 <- imodify(mat, ~ qunif(.x, umin[[.y]], umax[[.y]]))
mat2 <- mutate(mat2, across(everything(), exp), ID = row_number())
head(mat)
```

    . # A tibble: 6 × 5
    .    TVCL  TVVC  TVKA   TVQ  TVVP
    .   <dbl> <dbl> <dbl> <dbl> <dbl>
    . 1 0.5   0.5   0.5   0.5   0.5  
    . 2 0.75  0.25  0.75  0.25  0.75 
    . 3 0.25  0.75  0.25  0.75  0.25 
    . 4 0.375 0.375 0.625 0.125 0.875
    . 5 0.875 0.875 0.125 0.625 0.375
    . 6 0.625 0.125 0.375 0.375 0.125

# Run the analysis

## Simulation

``` r
out <- mrgsim_ei(mod, sunev(), mat2, output = "df")
out <- as.data.table(out)
```

## Calculate AUC

``` r
y <- out[, list(auc = auc_partial(time,CP)), by = "ID"][["auc"]]
```

## Indices

``` r
ind <- sobol_indices(Y = y, N = N, params = names(params), boot = TRUE, R = 1000, first = "jansen")
ind.dummy <- sobol_dummy(Y = y, N = N, params = names(params), boot = TRUE, R = 1000)
```

``` r
ind
```

    . 
    . First-order estimator: jansen | Total-order estimator: jansen 
    . 
    . Total number of model runs: 229376 
    . 
    . Sum of first order indices: 0.772696 
    .         original          bias    std.error       low.ci     high.ci
    .  1: 0.1800927710  3.713703e-04 0.0084472108  0.163165172 0.196277630
    .  2: 0.5169806384  4.025369e-04 0.0062918286  0.504246344 0.528909859
    .  3: 0.0722516683  4.572892e-04 0.0067438210  0.058576733 0.085012025
    .  4: 0.0029247659  4.538916e-04 0.0059300225 -0.009151756 0.014093505
    .  5: 0.0004461823  3.174973e-04 0.0055779353 -0.010803867 0.011061237
    .  6: 0.3828283356 -1.992348e-05 0.0057454181  0.371587447 0.394109072
    .  7: 0.7289651275 -5.446301e-04 0.0087529530  0.712354285 0.746665230
    .  8: 0.1162534021  6.669582e-05 0.0019249053  0.112413961 0.119959451
    .  9: 0.0119834278 -7.208326e-07 0.0003664219  0.011265975 0.012702322
    . 10: 0.0019679814 -1.318534e-06 0.0001110980  0.001751552 0.002187048
    .     sensitivity parameters
    .  1:          Si       TVCL
    .  2:          Si       TVVC
    .  3:          Si       TVKA
    .  4:          Si        TVQ
    .  5:          Si       TVVP
    .  6:          Ti       TVCL
    .  7:          Ti       TVVC
    .  8:          Ti       TVKA
    .  9:          Ti        TVQ
    . 10:          Ti       TVVP

## Visualize

``` r
plot(ind, dummy = ind.dummy) + ylim(0,1)
```

![](img/sensobolunnamed-chunk-14-1.png)<!-- -->
