function MZ =mf_mz_cs(S,X)

%This is stage 2 for Aligning Torque curve fitting 
%this contains coefficients related to Inclination Angle 

global V
for i=1:length(V)
    eval(strcat('global'," ",V{i,1}));
    eval(strcat(V{i,1},'=',num2str(V{i,2}),';'));
end
clear i

FZ     =  abs(X(:,1));
ALPHA  =  X(:,2)*pi/180;
KAPPA  =   X(:,3);
GAMMA  =  X(:,4)*pi/180;

GAMMAY = GAMMA .* LGAY; %31 (%48 lgay=lg
GAMMAZ = GAMMA .* LGAZ; %47 (%63 lgaz = lg
FZ0PR  = FZ0  .*  LFZO; %15,  NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ    = (FZ-FZ0PR) ./ FZ0PR; %14,  (%30)

SSZ1 = S(1);
SSZ2 = S(2);
SSZ3 = S(3);
SSZ4 = S(4);

ZETA1 = 1.00;
%C
%C -- ALIGNING TORQUE (PURE SIDE SLIP)
%C
SHY     = (PHY1+PHY2 .* DFZ) .* LHY + PHY3 .* GAMMAY; %38,  (%55)
SVY     = FZ .* ((PVY1+PVY2 .* DFZ) .* LVY+(PVY3+PVY4 .* DFZ) .* GAMMAY) .* LMUY; %39 (%56)
ALPHAY  = ALPHA+SHY;  %30 (%47)

SHT    = QHZ1+QHZ2 .* DFZ+(QHZ3+QHZ4 .* DFZ) .* GAMMAZ; %52 ( %68)

KX = FZ.*(PKX1 + PKX2.*DFZ).*exp(PKX3.*DFZ).* LKX;
KY      = PKY1 .* FZ0 .* sin(2.0 .* atan(FZ ./ (PKY2 .* FZ0 .* LFZO))) .* (1.0-PKY3 .* abs(GAMMAY)) .* LFZO .* LKY; %36 (%53)

ALPHAT = ALPHA+SHT;  %43 (%59)
ALPHAT = (ALPHAT.^2+((KX./KY).^2).*(KAPPA.^2)).*sign(ALPHAT);

% NOTE: PER SVEN, "EQUATION 45 IS WRONG DOCUMENTATION,
% THERE IT SHOULD BE SHF INSTEAD OF SHR"
SHF    = SHY+SVY ./ KY; %46 (%62)

ALPHAR = ALPHA+SHF; %45 (%61)
ALPHAR = (ALPHAR.^2+((KX./KY).^2).*(KAPPA.^2)).*sign(ALPHAR);

BT = (QBZ1+QBZ2 .* DFZ+QBZ3 .* DFZ.^2) .* (1.0+QBZ4 .* GAMMAZ+QBZ5 .* abs(GAMMAZ)) .* LKY ./ LMUY; %48 (%64)
CT = QCZ1; %49 (%65)
DT = FZ .* (QDZ1+QDZ2 .* DFZ) .* (1.0+QDZ3 .* GAMMAZ+QDZ4 .* GAMMAZ.^2) .* (R0 ./ FZ0) .* LTR; %50 (%66)
% NOTE: EQUATION FOR ET HAS CHANGED FROM PAC97 EQUATION; 2/PI TERM IS NEW.
ET = (QEZ1+QEZ2 .* DFZ+QEZ3 .* DFZ.^2) .* (1.0+(QEZ4+QEZ5 .* GAMMAZ) .* (2/pi) .* atan(BT .* CT .* ALPHAT)); %51 (%67)
CY      = PCY1 .* LCY;  %32 (%49)
MUY     = (PDY1+PDY2 .* DFZ) .* (1.0-PDY3 .* GAMMAY.^2) .* LMUY; %34 (%51)

DY      = MUY .* FZ; %33 (%50)
BY      = KY ./ (CY .* DY);  %37 (%54)
BR = QBZ9 .* LKY ./ LMUY+QBZ10 .* BY .* CY;  %53 (%69)
% NOTE: LRES MULTIPLIES EVERYTHING IN ORIG EQN; BELOW MATCHES DOCUMENTATION
DR = FZ .* ((QDZ6+QDZ7 .* DFZ) .* LRES+(QDZ8+QDZ9 .* DFZ) .* GAMMAZ) .* R0 .* LMUY; %54 (%70 LRES=LMR)

TRAIL = DT .* cos(CT .* atan(BT .* ALPHAT-ET .* (BT .* ALPHAT-atan(BT .* ALPHAT)))) .* cos(ALPHA); %42 (%58)
MZR = DR .* cos(atan(BR .* ALPHAR)) .* cos(ALPHA); %44 (%60)
EY      = (PEY1+PEY2 .* DFZ) .* (1.0-(PEY3+PEY4 .* GAMMAY) .* sign(ALPHAY)) .* LEY; %35 (%52)
% NOTE: LVY MULTIPLIES ONLY PVY1&2 IN DOCUMENTATION; ORIG VERSION MULT ALL TERMS
SVY     = FZ .* ((PVY1+PVY2 .* DFZ) .* LVY+(PVY3+PVY4 .* DFZ) .* GAMMAY) .* LMUY; %39 (%56)
FY0     = DY .* sin(CY .* atan(BY .* ALPHAY-EY .* (BY .* ALPHAY-atan(BY .* ALPHAY))))+SVY; %29 (%46)

%c
%c --Lateral Force (combined slip)
%c

BYK = (RBY1 + RBY4.* GAMMAY.* GAMMAY).* cos(atan(RBY2.*(ALPHA - RBY3)));

CYK = RCY1;

EYK = REY1 + REY2.* DFZ;

DVYK = MUY.*FZ.*(RVY1 + RVY2.* DFZ + RVY3.* GAMMAY).* cos(atan(RVY4.* ALPHA));

SHYK = RHY1 + RHY2.* DFZ;

SVYK = DVYK.*sin(RVY5.*atan(RVY6));

GYKO = cos(CYK.*atan(BYK.*SHYK - EYK.*(BYK.*SHYK - atan(BYK.*SHYK))));

KAPPA_S = KAPPA + SHYK;

GYK = cos(CYK.* atan(BYK.*KAPPA_S - EYK.*(BYK.*KAPPA_S - atan(BYK.*KAPPA_S))))./GYKO;

FY = GYK.* FY0 + SVYK;

S = R0.*(SSZ1 + SSZ2.*(FY/FZ0PR)+(SSZ3+SSZ4.*DFZ).*(GAMMAY.^2));



%c
%c -- longitudinal force (pure slip)
%c
SHX = (PHX1+PHX2 .* DFZ) .* LHX; %27,

KAPPAX = KAPPA + SHX; %19

GAMMAX = GAMMA * LGAY;

CX = PCX1 .* LCX; %21

MUX = (PDX1 + PDX2 .* DFZ) .* (1. -PDX3 .* GAMMAX.^2) .* LMUX; %23

DX = MUX .* FZ; %22

KX = FZ.*(PKX1 + PKX2.*DFZ).*exp(PKX3.*DFZ).* LKX;

BX = KX ./ (CX .* DX); %26
% NOTE, PER SVEN @TNO: "SIGN(ALPHAY)"IS CORRECT AS IN DOCUMENTATION & BELOW; IT'S NOT SUPPOSED TO BE "SIGN(GAMMAY)"

EX = (PEX1 + PEX2 .* DFZ +PEX3.*DFZ.*DFZ ) .* (1-PEX4 .*sign(KAPPAX)) ; %24
% NOTE: LVX MULTIPLIES ONLY PVX1&2 IN DOCUMENTATION; ORIG VERSION MULT ALL TERMS

SVX = FZ .* (PVX1 + PVX2 .* DFZ) .* LVX .* LMUX .* ZETA1; %28

FX0 = DX .* sin((CX .* atan(BX .* KAPPAX -EX .*(BX.*KAPPAX - atan(BX.*KAPPAX))))+SVX);

%--Longitudinal Force combined slip---

SHXA = RHX1;

CXA = RCX1;

EXA = REX1 + REX2.* DFZ;

BXA = (RBX1 + RBX3.* GAMMAX.* GAMMAX).* cos(atan(RBX2.* KAPPA));

ALPHA_S = sin(ALPHA) + SHXA;

GXA0 = cos(CXA.*atan(BXA.*SHXA - EXA.*(BXA.*SHXA - atan(BXA.*SHXA))));

GXA = cos(CXA.*atan(BXA.*ALPHA_S - EXA.*(BXA.* ALPHA_S - atan(BXA.*ALPHA_S))))./ GXA0;

FX = GXA .* FX0;



MZ= -TRAIL .* FY + MZR + S.*FX ;
