; Modified from appendix 4 of Jun Shen, Alison Boeckmann & Andrew Vick
; Modified for NONMEM 7.3
; Journal of Pharmacokinetics and Pharmacodynamics
; ISSN 1567-567X Volume 39 Number 3
; J Pharmacokinet Pharmacodyn (2012) 39:251-262 DOI 10.1007/s10928-012-9247-3

$PROBLEM Dose Superposition
$INPUT ID TIME CMT AMT DV ADDL II RATE
$DATA sumdosetf.csv IGNORE=@  ; identical to sall6.csv
$SUBROUTINES ADVAN6 TOL=10 TRANS1

$MODEL COMP=ABS ;(DEFDOSE)
       COMP=CENTRAL;(DEFOBS)
       COMP=PERI

$abbr declare dosetime(100),dose(100)
$abbr declare dowhile i
$abbr declare dowhile ndose

$PK 
CALLFL=-2
IF (NEWIND < 2) NDOSE=0

IF (AMT > 0 .and. cmt==1) THEN
 NDOSE=NDOSE+1
 dosetime(NDOSE)=TIME
 DOSE(NDOSE)=AMT
ENDIF

CL = THETA(1)*EXP(ETA(1))
V2 = THETA(2)*EXP(ETA(2))
Q  = THETA(3)
V3 = THETA(4)

K=CL/V2
S2 =V2
K23=Q/V2
K32=Q/V3

;--------ABSORPTION MODEL-----------
F1=0
MTT  = THETA(5)*EXP(ETA(4))
NN   = THETA(6)*EXP(ETA(5))
BIO  = THETA(7)*EXP(ETA(3))
KTR  = (NN+1)/MTT

NFAC=SQRT(2*3.1416)*NN**(NN+0.5)*(EXP(-NN))*(1+1/(12*NN))
KINPT=BIO*KTR**(NN+1)/NFAC

$DES
INPT=0
I=1
DOWHILE (I<=NDOSE)
IPT=0
IF (T>=dosetime(I)) IPT=DOSE(I)*(T-dosetime(I))**NN*EXP(-KTR*(T-dosetime(I)))
INPT=INPT+IPT
I=I+1
ENDDO

 DADT(1)=KINPT*INPT-KTR*A(1)
 DADT(2)=KTR*A(1)-K23*A(2)-K*A(2)+K32*A(3)
 DADT(3)=K23*A(2)-K32*A(3)

$ERROR

Y = F*EXP(EPS(1))

$THETA 
4 ;CL
3 ;V2
1 ;Q
2 ;V3
10 ;MTT
5 ;NN
0.85 ;BIO

$OMEGA 0 FIX 0 FIX 0 FIX 0 FIX 0 FIX 
$SIGMA 0 FIX
$SIMULATION (123456) ONLY
;$ESTIMATION METHOD=0
$TABLE  ID TIME CMT CL V2 Q V3 MTT NN INPT ONEHEADER NOPRINT FORMAT=s1PE18.11 FILE=sumdosetn.txt