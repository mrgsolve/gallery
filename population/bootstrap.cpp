$PROB
# Example population PK model

$PARAM  WT=70

$THETA
log(1)  log(24) log(0.5)

$PKMODEL cmt="GUT CENT", depot=TRUE

$SET end=240, delta=0.5

$MAIN
double CL = exp(THETA1  + 0.75*log(WT/70) + ECL);
double V  = exp(THETA2  +      log(WT/70) + EV );
double KA = exp(THETA3                    + EKA);

$OMEGA @labels ECL EV EKA
0.3 0.1 0.5

$SIGMA 0 

$TABLE
capture IPRED = CENT/V;
capture DV = IPRED*exp(EPS(1));

$CAPTURE CL V ECL

  
