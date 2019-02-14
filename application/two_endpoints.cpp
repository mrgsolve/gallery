[ SET ] end = 72, delta = 0.25, add = 0.05

[ PARAM ] CLnr = 0.97, V = 22.3, KA = 1.9, CLr = 0.2, dvtype = 0

[ CMT ] GUT CENT URINE

[ ODE ] 
dxdt_GUT   = -KA*GUT;
dxdt_CENT  =  KA*GUT - CLnr*(CENT/V) - CLr*(CENT/V);
dxdt_URINE =  CLr*(CENT/V); 

[ SIGMA ] 0 0

[ TABLE ] 
capture CP = (CENT/V)*exp(EPS(1));
capture UR = URINE*exp(EPS(2));
capture DV = dvtype ==2 ? UR : CP;
