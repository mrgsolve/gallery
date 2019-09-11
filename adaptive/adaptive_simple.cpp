[ prob ]
Simple code to adjust doses.

[ param ]
CL = 1, V = 20, KA = 1.1, F1 = 1, F1adjust = 0.5, condition = 100

[ pkmodel ] cmt = "DEPOT CENT", depot = TRUE

[ preamble ] 
bool condition_met = false;

[ main ]
if(NEWIND <=1) condition_met = false;

F_DEPOT = F1;

if(condition_met) F_DEPOT = F1 * F1adjust;

[ table ] 
capture CP = CENT/V;

condition_met = CP > condition || condition_met;

[ capture ] F_DEPOT condition_met

