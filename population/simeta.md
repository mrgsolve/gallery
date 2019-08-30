Resimulate random effect variates on demand
================

# simeta (and simeps)

``` r
library(mrgsolve)
library(dplyr)
```

## `simeta` example

  - In this example, we want to simulate a patient-specific baseline
    response that is between 80 and 120.
  - In the code, we start a loop that calls `simeta` with no arguments
    until the baseline is between the specified bounds
  - For this example, we only calculate `BASE` when `NEWIND <=1` … or
    whenever we are working on the first record of an individual. This
    ensures that we don’t re-simulate `BASE` at every simulation record.
  - We have also implemented a counter (`i`) that ensures we only try to
    resimulate 100 times. If a value for `BASE` cannot be generated in
    100 tries, we give up.  
  - You probably won’t need to implement `FLAG` for your problem. I am
    only using `FLAG` here so we can selectively call the `simeta` code
    to demonstrate how it is working.

<!-- end list -->

``` r
code <- '
$PARAM TVBASE = 100, FLAG = 0

$CMT RESPONSE

$MAIN 

if(NEWIND <=1) {

  capture BASE = TVBASE*exp(EBASE);

  int i = 0;

  if(FLAG > 0) {
    while((BASE < 80) || (BASE > 120)) {
      if(++i > 100) {
        mrgsolve::report("There was a problem simulating BASE");
      }
      simeta();
      BASE = TVBASE*exp(EBASE);
    }
  }
  
  RESPONSE_0 = BASE;
}


$OMEGA @labels EBASE
1

$CAPTURE EBASE
'
```

``` r
mod <- mcode("simeta", code)
```

### Simulate without `simeta`

``` r
system.time({
  out <- mod %>% mrgsim(nid=100, end=-1)
  sum <- summary(out)
})
```

    .    user  system elapsed 
    .   0.021   0.001   0.022

``` r
print(sum)
```

    .        ID              time      RESPONSE           EBASE         
    .  Min.   :  1.00   Min.   :0   Min.   :  3.681   Min.   :-3.30205  
    .  1st Qu.: 25.75   1st Qu.:0   1st Qu.: 47.744   1st Qu.:-0.73971  
    .  Median : 50.50   Median :0   Median :123.320   Median : 0.20961  
    .  Mean   : 50.50   Mean   :0   Mean   :155.289   Mean   :-0.02367  
    .  3rd Qu.: 75.25   3rd Qu.:0   3rd Qu.:179.980   3rd Qu.: 0.58768  
    .  Max.   :100.00   Max.   :0   Max.   :846.886   Max.   : 2.13640  
    .       BASE        
    .  Min.   :  3.681  
    .  1st Qu.: 47.744  
    .  Median :123.320  
    .  Mean   :155.289  
    .  3rd Qu.:179.980  
    .  Max.   :846.886

When we simulate with `FLAG=0`, the `simeta` code **isn’t** called and
we `BASE` values all over the place.

### Simulate with `simeta`

``` r
system.time({
  out <- mod %>% mrgsim(nid=100, end=-1, param=list(FLAG=1))
  sum <- summary(out)
})
```

    .    user  system elapsed 
    .   0.004   0.000   0.004

``` r
print(sum)
```

    .        ID              time      RESPONSE          EBASE         
    .  Min.   :  1.00   Min.   :0   Min.   : 80.45   Min.   :-0.21758  
    .  1st Qu.: 25.75   1st Qu.:0   1st Qu.: 86.76   1st Qu.:-0.14205  
    .  Median : 50.50   Median :0   Median : 98.83   Median :-0.01179  
    .  Mean   : 50.50   Mean   :0   Mean   : 98.46   Mean   :-0.02222  
    .  3rd Qu.: 75.25   3rd Qu.:0   3rd Qu.:108.18   3rd Qu.: 0.07862  
    .  Max.   :100.00   Max.   :0   Max.   :119.91   Max.   : 0.18160  
    .       BASE       
    .  Min.   : 80.45  
    .  1st Qu.: 86.76  
    .  Median : 98.83  
    .  Mean   : 98.46  
    .  3rd Qu.:108.18  
    .  Max.   :119.91

When we simulate with `FLAG=1`, the `simeta` code **is** called and we
`BASE` values within the specified bounds.

## `simeps` example

  - In this example, we want to re-simulate the residual error variate
    to make sure we have a concentration that is positive.
  - We set up a loop that looks like the `simeta` example, but we work
    in `$TABLE` this time because we are calculating `CP`.

<!-- end list -->

``` r
code <- '
$PARAM CL = 1, V = 20, FLAG=0

$SIGMA 50

$PKMODEL cmt="CENT"

$TABLE
capture CP = CENT/V + EPS(1);

int i = 0;
while(CP < 0 && FLAG > 0) {
  if(++i > 100) {
    mrgsolve::report("Problem simulating positive CP");
  }
  simeps();
  CP = CENT/V + EPS(1);
}

'
```

``` r
mod <- mcode("simeps", code)
```

### Simulate without `simeps`

``` r
system.time({
  out <- mod %>% ev(amt=100) %>% mrgsim(end=48)
  sum <- summary(out)
})
```

    .    user  system elapsed 
    .   0.011   0.000   0.011

``` r
print(sum)
```

    .        ID         time            CENT              CP         
    .  Min.   :1   Min.   : 0.00   Min.   :  0.00   Min.   :-13.935  
    .  1st Qu.:1   1st Qu.:11.25   1st Qu.: 15.93   1st Qu.: -3.445  
    .  Median :1   Median :23.50   Median : 29.38   Median :  1.410  
    .  Mean   :1   Mean   :23.52   Mean   : 37.47   Mean   :  1.846  
    .  3rd Qu.:1   3rd Qu.:35.75   3rd Qu.: 54.21   3rd Qu.:  5.705  
    .  Max.   :1   Max.   :48.00   Max.   :100.00   Max.   : 21.803

**Negative** concentrations are simulated when we **don’t** call the
`simeps` loop.

### Simulate with `simeps`

``` r
system.time({
  out <- mod %>% ev(amt=100) %>% mrgsim(end=48, param=list(FLAG=1))
  sum <- summary(out)
})
```

    .    user  system elapsed 
    .   0.032   0.000   0.033

``` r
print(sum)
```

    .        ID         time            CENT              CP        
    .  Min.   :1   Min.   : 0.00   Min.   :  0.00   Min.   : 0.375  
    .  1st Qu.:1   1st Qu.:11.25   1st Qu.: 15.93   1st Qu.: 3.226  
    .  Median :1   Median :23.50   Median : 29.38   Median : 5.950  
    .  Mean   :1   Mean   :23.52   Mean   : 37.47   Mean   : 6.430  
    .  3rd Qu.:1   3rd Qu.:35.75   3rd Qu.: 54.21   3rd Qu.: 8.263  
    .  Max.   :1   Max.   :48.00   Max.   :100.00   Max.   :21.632

Better … all concentrations are positive.
