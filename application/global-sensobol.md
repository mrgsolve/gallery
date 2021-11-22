Sobol sensitivity analysis using sensobol
================
Kyle Baron
2021-11-22 16:27:20

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

For this example, we will just take uniform samples based on the the
current value of the parameter.

``` r
params <- unlist(as.list(param(mod))[c("TVCL", "TVVC", "TVKA", "TVQ", "TVVP")])
umin <- params / 5
umax <- params * 5
umin
```

    .    TVCL    TVVC    TVKA     TVQ    TVVP 
    .  10.360 406.000   0.039   1.444 116.600

After setting the minimum and maximum for each parameter, get the value
of the parameter by looking at the quantiles of the (log uniform)
distribution

``` r
mat2 <- imodify(mat, ~ qunif(.x, umin[[.y]], umax[[.y]]))
mat2 <- mutate(mat2, ID = row_number())
head(mat2)
```

    . # A tibble: 6 × 6
    .    TVCL  TVVC  TVKA   TVQ  TVVP    ID
    .   <dbl> <dbl> <dbl> <dbl> <dbl> <int>
    . 1 135.   5278 0.507 18.8  1516.     1
    . 2 197.   2842 0.741 10.1  2215.     2
    . 3  72.5  7714 0.273 27.4   816.     3
    . 4 104.   4060 0.624  5.78 2565.     4
    . 5 228.   8932 0.156 23.1  1166      5
    . 6 166.   1624 0.39  14.4   466.     6

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
    . Sum of first order indices: 0.7380731 
    .          original          bias    std.error       low.ci     high.ci
    .  1:  0.1386450195 -5.443701e-05 0.0134605344  0.112317294 0.165081619
    .  2:  0.5702287103 -1.415956e-04 0.0102223164  0.550334934 0.590405678
    .  3:  0.0302745338  2.545441e-04 0.0061872971  0.017893110 0.042146869
    .  4:  0.0000557737  1.652504e-04 0.0063841189 -0.012622120 0.012403167
    .  5: -0.0011309726 -4.711711e-05 0.0054836856 -0.011831682 0.009663971
    .  6:  0.3885564747  1.018732e-04 0.0100190499  0.368817625 0.408091579
    .  7:  0.8235786500  6.278763e-04 0.0138431873  0.795818625 0.850082922
    .  8:  0.0429891387  2.074858e-05 0.0012776991  0.040464146 0.045472634
    .  9:  0.0084979002 -6.232607e-06 0.0005036811  0.007516936 0.009491330
    . 10:  0.0013766977 -2.950487e-06 0.0001292514  0.001126320 0.001632976
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

Plot outputs versus inputs

``` r
plot_scatter(N = 2000, data = mat2, Y = y, params = names(mat))
```

![](img/sensobolunnamed-chunk-15-1.png)<!-- -->

``` r
plot_multiscatter(N = 2000, data = mat2, Y = y, params = names(mat))
```

![](img/sensobolunnamed-chunk-16-1.png)<!-- -->
