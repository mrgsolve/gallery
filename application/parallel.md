mrgsolve in Parallel
================
Kyle Baron
2018-06-25 23:15:53

-   [About](#about)
-   [An example model](#an-example-model)
-   [The `future.apply` package](#the-future.apply-package)
    -   [Choose a plan](#choose-a-plan)
    -   [Simulate with `future_lapply`](#simulate-with-future_lapply)
-   [Compare methods](#compare-methods)
    -   [`future_lapply`](#future_lapply)
    -   [`lapply`](#lapply)
    -   [`mclapply`](#mclapply)

About
=====

This vignette looks at options for parallelizing simulations with mrgsolve in a platform-independent way. We utilize the `future.apply` package (available on CRAN) to do this.

Your mileage may vary in terms of speedup factor. It is highly dependent on the problem you have. Also, with any method there is some overhead that needs to be taken into consideration when planning the simulations. It is very possible that your parallelized setup takes **longer** with the non-parallel setup.

An example model
================

``` r
library(dplyr)
library(mrgsolve)
mod <- mread("pk1", modlib())
```

The `future.apply` package
==========================

``` r
library(future.apply)
```

Choose a plan
-------------

``` r
plan("multicore")
```

Simulate with `future_lapply`
-----------------------------

Works pretty much like `lapply`

``` r
out <- future_lapply(1:10, function(i) {
  mod %>% 
    ev(amt = 100) %>%
    mrgsim() %>% 
    mutate(i = i)
}) %>% bind_rows
```

``` r
head(out)
```

    . # A tibble: 6 x 6
    .      ID  time     EV  CENT    CP     i
    .   <dbl> <dbl>  <dbl> <dbl> <dbl> <int>
    . 1     1     0   0      0    0        1
    . 2     1     0 100      0    0        1
    . 3     1     1  36.8   61.4  3.07     1
    . 4     1     2  13.5   81.0  4.05     1
    . 5     1     3   4.98  85.4  4.27     1
    . 6     1     4   1.83  84.3  4.21     1

Compare methods
===============

`future_lapply`
---------------

``` r
system.time({
  out <- future_lapply(1:2000, function(i) {
    mod %>% 
      ev(amt = 100, ii = 24, addl = 27) %>%
      mrgsim(end = 28*24, nid = 20) %>% 
      mutate(i = i)
  }) %>% bind_rows
})
```

    .    user  system elapsed 
    .  36.781   6.721  11.518

`lapply`
--------

``` r
system.time({
  out <- lapply(1:2000, function(i) {
    mod %>% 
      ev(amt = 100, ii = 24, addl = 27) %>%
      mrgsim(end = 28*24, nid = 20) %>% 
      mutate(i = i)
  }) %>% bind_rows
})
```

    .    user  system elapsed 
    .  19.613   1.246  21.278

`mclapply`
----------

``` r
system.time({
  out <- parallel::mclapply(1:2000, function(i) {
    mod %>% 
      ev(amt = 100, ii = 24, addl = 27) %>%
      mrgsim(end = 28*24, nid = 20) %>% 
      mutate(i = i)
  }) %>% bind_rows
})
```

    .    user  system elapsed 
    .  21.460   6.112  17.194
