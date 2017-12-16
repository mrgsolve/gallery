Resimulate random effect variates on demand
================

simeta (and simeps)
===================

``` r
library(mrgsolve)
library(dplyr)
```

`simeta` example
----------------

-   In this example, we want to simulate a patient-specific baseline response that is between 80 and 120.
-   In the code, we start a loop that calls `simeta` with no arguments until the baseline is between the specified bounds
-   For this example, we only calculate `BASE` when `NEWIND <=1` ... or whenever we are working on the first record of an individual. This ensures that we don't re-simulate `BASE` at every simulation record.
-   We have also implemented a counter (`i`) that ensures we only try to resimulate 100 times. If a value for `BASE` cannot be generated in 100 tries, we give up.
-   You probably won't need to implement `FLAG` for your problem. I am only using `FLAG` here so we can selectively call the `simeta` code to demonstrate how it is working.

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
        report("There was a problem simulating BASE");
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
    .   0.010   0.001   0.014

``` r
print(sum)
```

    .        ID              time      RESPONSE          EBASE         
    .  Min.   :  1.00   Min.   :0   Min.   : 10.96   Min.   :-2.21128  
    .  1st Qu.: 25.75   1st Qu.:0   1st Qu.: 55.25   1st Qu.:-0.59340  
    .  Median : 50.50   Median :0   Median :104.56   Median : 0.04464  
    .  Mean   : 50.50   Mean   :0   Mean   :168.07   Mean   : 0.08482  
    .  3rd Qu.: 75.25   3rd Qu.:0   3rd Qu.:242.55   3rd Qu.: 0.88605  
    .  Max.   :100.00   Max.   :0   Max.   :914.65   Max.   : 2.21337  
    .       BASE       
    .  Min.   : 10.96  
    .  1st Qu.: 55.25  
    .  Median :104.56  
    .  Mean   :168.07  
    .  3rd Qu.:242.55  
    .  Max.   :914.65

When we simulate with `FLAG=0`, the `simeta` code **isn't** called and we `BASE` values all over the place.

### Simulate with `simeta`

``` r
system.time({
  out <- mod %>% mrgsim(nid=100, end=-1, param=list(FLAG=1))
  sum <- summary(out)
})
```

    .    user  system elapsed 
    .   0.102   0.005   0.108

``` r
print(sum)
```

    .        ID              time      RESPONSE          EBASE         
    .  Min.   :  1.00   Min.   :0   Min.   : 80.11   Min.   :-0.22175  
    .  1st Qu.: 25.75   1st Qu.:0   1st Qu.: 88.56   1st Qu.:-0.12147  
    .  Median : 50.50   Median :0   Median : 98.66   Median :-0.01345  
    .  Mean   : 50.50   Mean   :0   Mean   : 98.92   Mean   :-0.01713  
    .  3rd Qu.: 75.25   3rd Qu.:0   3rd Qu.:108.45   3rd Qu.: 0.08115  
    .  Max.   :100.00   Max.   :0   Max.   :119.78   Max.   : 0.18047  
    .       BASE       
    .  Min.   : 80.11  
    .  1st Qu.: 88.56  
    .  Median : 98.66  
    .  Mean   : 98.92  
    .  3rd Qu.:108.45  
    .  Max.   :119.78

When we simulate with `FLAG=1`, the `simeta` code **is** called and we `BASE` values within the specified bounds.

`simeps` example
----------------

-   In this example, we want to re-simulate the residual error variate to make sure we have a concentration that is positive.
-   We set up a loop that looks like the `simeta` example, but we work in `$TABLE` this time because we are calculating `CP`.

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
    report("Problem simulating positive CP");
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
    .   0.008   0.001   0.012

``` r
print(sum)
```

    .        ID         time            CENT              CP         
    .  Min.   :1   Min.   : 0.00   Min.   :  0.00   Min.   :-21.477  
    .  1st Qu.:1   1st Qu.:11.25   1st Qu.: 15.93   1st Qu.: -5.464  
    .  Median :1   Median :23.50   Median : 29.38   Median :  2.901  
    .  Mean   :1   Mean   :23.52   Mean   : 37.47   Mean   :  0.224  
    .  3rd Qu.:1   3rd Qu.:35.75   3rd Qu.: 54.21   3rd Qu.:  6.046  
    .  Max.   :1   Max.   :48.00   Max.   :100.00   Max.   : 18.789

**Negative** concentrations are simulated when we **don't** call the `simeps` loop.

### Simulate with `simeps`

``` r
system.time({
  out <- mod %>% ev(amt=100) %>% mrgsim(end=48, param=list(FLAG=1))
  sum <- summary(out)
})
```

    .    user  system elapsed 
    .   0.048   0.003   0.052

``` r
print(sum)
```

    .        ID         time            CENT              CP         
    .  Min.   :1   Min.   : 0.00   Min.   :  0.00   Min.   : 0.2136  
    .  1st Qu.:1   1st Qu.:11.25   1st Qu.: 15.93   1st Qu.: 2.1379  
    .  Median :1   Median :23.50   Median : 29.38   Median : 4.7224  
    .  Mean   :1   Mean   :23.52   Mean   : 37.47   Mean   : 5.4664  
    .  3rd Qu.:1   3rd Qu.:35.75   3rd Qu.: 54.21   3rd Qu.: 8.7779  
    .  Max.   :1   Max.   :48.00   Max.   :100.00   Max.   :16.5729

Better ... all concentrations are positive.
