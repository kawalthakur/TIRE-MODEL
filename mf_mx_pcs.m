function MX =mf_mx_pcs(C,X)

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
FZ0PR  = FZ0  .*  LFZO; %15,  NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ    = (FZ-FZ0PR) ./ FZ0PR; %14,  (%30)

QSX1    = C(1);
QSX2    = C(2);
QSX3    = C(3);
QSX4    = C(4);
QSX5    = C(5);
QSX6    = C(6);
QSX7    = C(7);
QSX8    = C(8);
QSX9    = C(9);
QSX10   = C(10);
QSX11   = C(11);

%c
%c -- lateral force (pure side slip)
%c
SHY     = (PHY1+PHY2 .* DFZ) .* LHY + PHY3 .* GAMMAY; %38,  (%55)
ALPHAY  = ALPHA+SHY;  %30 (%47)
CY      = PCY1 .* LCY;  %32 (%49)
MUY     = (PDY1+PDY2 .* DFZ) .* (1.0-PDY3 .* GAMMAY.^2) .* LMUY; %34 (%51)
DY      = MUY .* FZ; %33 (%50)
KY      = PKY1 .* FZ0 .* sin(2.0 .* atan(FZ ./ (PKY2 .* FZ0 .* LFZO))) .* (1.0-PKY3 .* abs(GAMMAY)) .* LFZO .* LKY; %36 (%53)
BY      = KY ./ (CY .* DY);  %37 (%54)
% NOTE, PER SVEN @TNO: "SIGN(ALPHAY)"IS CORRECT AS IN DOCUMENTATION & BELOW; IT'S NOT SUPPOSED TO BE "SIGN(GAMMAY)"
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

DVYK = MUY.*FZ.*(RVY1 + RVY2.* DFZ + RVY3.* GAMMA).* cos(atan(RVY4.* ALPHA));

SHYK = RHY1 + RHY2.* DFZ;

SVYK = DVYK.*sin(RVY5.*atan(RVY6));

GYKO = cos(CYK.*atan(BYK.*SHYK - EYK.*(BYK.*SHYK - atan(BYK.*SHYK))));

KAPPA_S = KAPPA + SHYK;

GYK = cos(CYK.* atan(BYK.*KAPPA_S - EYK.*(BYK.*KAPPA_S - atan(BYK.*KAPPA_S))))./GYKO;

FY = GYK.* FY0 + SVYK;




MX = FZ.*R0.*((QSX1-QSX2.*GAMMA+QSX3.*FY./FZ0)+(QSX3.*FZ./FZ0)...
    +(QSX4.*cos(QSX5.*atan((QSX6.*FZ/FZ0).^2)).*sin((QSX7.*GAMMA)+QSX8*atan(QSX9.*FY./FZ0)))...
    +(QSX10.*atan(QSX11.*FZ./FZ0).*GAMMA)).*LMUX;      %84  SWIFT DOES NOT IMPLEMENT MX

