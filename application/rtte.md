Simulating large time-to-event trials with mrgsolve
================

# Introduction

This vignette shows implementation of MTIME approach to simulating large
time-to-event trial in NONMEM. We show the Repeated Time-to-Event model.

## Reference

**Simulating large time-to-event trials in NONMEM**

Joakim Nyberg, Kristin E. Karlsson, Siv Jönsson, Ulrika S.H. Simonsson,
Mats O. Karlsson and Andrew C. Hooker

PAGE 23 (2014) Abstr 3166 <https://www.page-meeting.org/?abstract=3166>

``` r
library(tidyverse)
theme_set(theme_bw())
library(mrgsolve)
```

# Input data

``` r
data <- data.frame(
  TIME = c(0, 288), 
  CMT = 0, 
  AMT = 0, 
  EVID = 0
)
ids <- seq(1, 300)
data <- crossing(ID = ids, data)

head(data)
```

    . # A tibble: 6 × 5
    .      ID  TIME   CMT   AMT  EVID
    .   <int> <dbl> <dbl> <dbl> <dbl>
    . 1     1     0     0     0     0
    . 2     1   288     0     0     0
    . 3     2     0     0     0     0
    . 4     2   288     0     0     0
    . 5     3     0     0     0     0
    . 6     3   288     0     0     0

# RTTE Model

``` r
mod <- mread("rtte.txt", end = -1)
```

    . Building rtte_txt ... done.

``` r
set.seed(98765)
out <- mrgsim(mod, data)

head(out)
```

    .   ID TIME       A1
    . 1  1    0 0.000000
    . 2  1  288 2.155004
    . 3  2    0 0.000000
    . 4  2  288 1.395972
    . 5  3    0 0.000000
    . 6  3  288 1.003641

# Check event output

``` r
res <- as.list(mod@envir)
x <- matrix(res$data, ncol = length(res$names), byrow=TRUE)
x <- as.data.frame(x)
names(x) <- res$names

head(x)
```

    .   ID DV       TIME         R        BASE
    . 1  1  1 111.989477 0.4327388 0.007482653
    . 2  1  1 175.056369 0.6240704 0.007482653
    . 3  2  1 104.784654 0.6018752 0.004847126
    . 4  4  1   9.111563 0.9376066 0.007090281
    . 5  4  1 238.529145 0.1966469 0.007090281
    . 6  4  1 241.393985 0.9804362 0.007090281

## Calculate analytic Event Time

``` r
x <- mutate(
  x, 
  AT = -log(R)/BASE, 
  DIFF = TIME - AT
)

head(x)
```

    ##   ID DV       TIME         R        BASE         AT         DIFF
    ## 1  1  1 111.989477 0.4327388 0.007482653 111.941713   0.04776402
    ## 2  1  1 175.056369 0.6240704 0.007482653  63.011350 112.04501842
    ## 3  2  1 104.784654 0.6018752 0.004847126 104.743544   0.04111040
    ## 4  4  1   9.111563 0.9376066 0.007090281   9.086349   0.02521333
    ## 5  4  1 238.529145 0.1966469 0.007090281 229.376715   9.15242952
    ## 6  4  1 241.393985 0.9804362 0.007090281   2.786591 238.60739373

## Difference between MTIME and Analytic Time

Take time to first event

``` r
xx <- distinct(x, ID, .keep_all = TRUE)
ggplot(xx, aes(x = DIFF)) + 
  geom_histogram(color = "white")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](rtte_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Model code

``` r
[ theta ] 0.0058 0.1

[ omega ] 0.09

[ set ] rtol = 1e-9, atol = 1e-9

[ plugin ]  autodec,  Rcpp, mrgx

[ cmt ] A1

[ global ] 
std::vector<double> COM(9);
std::vector<double> out;

[ pk ] 
BASE = THETA(1) * exp(ETA(1));

if(NEWIND==0) {
  COM[6] = 1;
  out.clear();
}

if(NEWIND==1) {
  ICOUNT = COM[6] + 1;
  COM[6] = ICOUNT;
}

if(NEWIND != 2) {
  RR = R::runif(0,1);
  COM[4] = 288; 
  COM[3] = -1; 
  COM[2] = RR; 
  COM[1] = -1;
  COM[7] = 0; 
  COM[8] = 0;
}

id = ID;

[ des ] 
dxdt_A1 = BASE;

SUR = exp(-(A1-COM[8]));
XR = 0;

if(COM[2] > SUR) {
  COM[1] = SOLVERTIME;
  COM[3] = SUR; 
  COM[8] = A1;
  COM[7] = COM[7] + 1;
  RR = R::runif(0,1);
  TMP = COM[2];
  COM[2] = RR;
  MYDV = 1;
  if(SOLVERTIME < COM[4]) {
    out.push_back(id); 
    out.push_back(MYDV); 
    out.push_back(COM[1]);
    out.push_back(TMP);
    out.push_back(BASE);
  }
}

[ error ] 

if(TIME <= 288) {
  mt = self.mtime(TIME + THETA(2));
}

if(self.rown +1 == self.nrow) {
  std::vector<std::string> names; 
  names.push_back("ID");
  names.push_back("DV"); 
  names.push_back("TIME"); 
  names.push_back("R"); 
  names.push_back("BASE");
  mrgx::get_envir(self).assign("data", out);
  mrgx::get_envir(self).assign("names", names);
}
```
