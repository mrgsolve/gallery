Sobol sensitivity analysis using sensobol
================
Kyle Baron
2021-11-22 16:05:53

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

We are just looking at a single dose for now.

``` r
sunev <- function(amt = 50,...) ev(amt = amt, ...)
```

# Generate samples

**Number**

``` r
N <- 2 ^ 15
```

**Generate**

``` r
mat <- sobol_matrices(N = N, params = c("TVCL", "TVVC", "TVKA", "TVQ", "TVVP"))
mat <- as_tibble(mat)
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

**Transform and groom**

For this example, we will just take uniform samples based on the
logarithm of the current value of the parameter.

``` r
params <- unlist(as.list(param(mod))[c("TVCL", "TVVC", "TVKA", "TVQ", "TVVP")])
umin <- log(params / 5)
umax <- log(params * 5)
umin
```

    .      TVCL      TVVC      TVKA       TVQ      TVVP 
    .  2.337952  6.006353 -3.244194  0.367417  4.758749

After setting the minimum and maximum for each parameter, get the value
of the parameter by looking at the quantiles of the (log uniform)
distribution

``` r
mat2 <- imodify(mat, ~ qunif(.x, umin[[.y]], umax[[.y]]))
mat2 <- mutate(mat2, across(everything(), exp), ID = row_number())
head(mat2)
```

    . # A tibble: 6 × 6
    .    TVCL  TVVC   TVKA   TVQ  TVVP    ID
    .   <dbl> <dbl>  <dbl> <dbl> <dbl> <int>
    . 1  51.8 2030  0.195   7.22  583      1
    . 2 116.   908. 0.436   3.23 1304.     2
    . 3  23.2 4539. 0.0872 16.1   261.     3
    . 4  34.6 1358. 0.292   2.16 1949.     4
    . 5 173.  6788. 0.0583 10.8   390.     5
    . 6  77.5  607. 0.130   4.83  174.     6

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
    .  1: 0.1800927710  1.238243e-05 0.0081723319  0.164062912 0.196097865
    .  2: 0.5169806384 -7.072754e-06 0.0064265026  0.504391998 0.529583425
    .  3: 0.0722516683  3.017526e-05 0.0065358709  0.059411422 0.085031565
    .  4: 0.0029247659  3.563857e-05 0.0057386288 -0.008358378 0.014136633
    .  5: 0.0004461823  7.463826e-05 0.0053585629 -0.010131046 0.010874134
    .  6: 0.3828283356  3.441830e-04 0.0059880652  0.370747760 0.394220545
    .  7: 0.7289651275  1.251542e-04 0.0082027993  0.712762782 0.744917165
    .  8: 0.1162534021  1.782163e-05 0.0019111638  0.112489768 0.119981393
    .  9: 0.0119834278  1.625625e-05 0.0003719938  0.011238077 0.012696266
    . 10: 0.0019679814  2.644654e-06 0.0001137576  0.001742376 0.002188298
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
