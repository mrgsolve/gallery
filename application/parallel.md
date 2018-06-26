mrgsolve in Parallel
================
Kyle Baron
2018-06-26 09:31:40

-   [About](#about)
-   [An example model](#an-example-model)
-   [The `future.apply` package](#the-future.apply-package)
    -   [Simulate with `future_lapply`](#simulate-with-future_lapply)
-   [Compare methods](#compare-methods)
    -   [`future_lapply`](#future_lapply)
    -   [`lapply`](#lapply)
    -   [`mclapply`](#mclapply)
-   [Parallelize within data set](#parallelize-within-data-set)

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

Simulate with `future_lapply`
-----------------------------

Works pretty much like `lapply`

``` r
plan("multiprocess")
```

Note: with `plan(multiprocess)`, you have to load the model shared object into the process. See `?laodso`.

``` r
e <- ev(amt = 100)
end <- 24

out <- future_lapply(1:10, function(i) {
  
  loadso(mod) ## NOTE
  
  mod %>% 
    ev(e) %>%
    mrgsim(end = end) %>% 
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

On macos or unix systems, you can use:

``` r
plan("multicore")
```

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
plan("multicore")

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
    .  36.205   6.555   9.204

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
    .  17.262   1.019  18.371

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
    .  19.887   4.320  13.706

Parallelize within data set
===========================

In this example, let's simulate 3k subjects at each of 8 doses. We'll split the data set on the dose and simulate each dose separately and then bind back together in a single data set. This is probably the quickest way to get it done. But we really need to work to see the speedup from parallelizing.

``` r
data <- expand.ev(
  ID = seq(2000), 
  amt = c(1,3,10,30,100,300,1000,3000),
  ii = 24, addl = 27
) 
count(data,amt)
```

    . # A tibble: 8 x 2
    .     amt     n
    .   <dbl> <int>
    . 1     1  2000
    . 2     3  2000
    . 3    10  2000
    . 4    30  2000
    . 5   100  2000
    . 6   300  2000
    . 7  1000  2000
    . 8  3000  2000

``` r
data_split <- split(data, data$amt)

system.time({
  out <- future_lapply(data_split, function(chunk) {
    mod %>% mrgsim_d(chunk, end = 24*27) %>% as_data_frame  
  }) %>% bind_rows()
})
```

    .    user  system elapsed 
    .   6.026   2.435   2.320

``` r
dim(out)
```

    . [1] 10400000        5

``` r
system.time({
  out <- lapply(data_split, function(chunk) {
    mod %>% mrgsim_d(chunk, end = 24*27) %>% as_data_frame  
  }) %>% bind_rows()
})
```

    .    user  system elapsed 
    .   3.326   0.460   3.805
