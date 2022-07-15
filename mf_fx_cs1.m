function FX = mf_fx_cs1(R,X)

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

%Coefficient for combined slip
RHX1 = R(1);
RCX1 = R(2);
REX1 = R(3);
REX2 = R(4);
RBX1 = R(5);
RBX2 = R(6);
%RBX3 = R(7);

% FZ0 = A(16) .* LFZO; %15, NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE

DFZ = (FZ-(FZ0*LFZO)) ./ (FZ0*LFZO) ; %14, (%30)

ZETA1 = 1;
LEX = 1;
LXA = 1;
%c
%c -- longitudinal force (pure slip)
%c
SHX = (PHX1+PHX2 .* DFZ) .* LHX; %27,

KAPPAX = KAPPA + SHX; %19

GAMMAX = GAMMA * LGAY;

CX = PCX1 .* LCX; %21

MUX = (PDX1 + PDX2 .* DFZ) .* (1. -PDX3 .* GAMMAX.^2) .* LMUX; %23

DX = MUX .* FZ .* ZETA1; %22

KX = FZ.*(PKX1 + PKX2.*DFZ).*exp(PKX3.*DFZ).* LKX;

BX = KX ./ (CX .* DX); %26
% NOTE, PER SVEN @TNO: "SIGN(ALPHAY)"IS CORRECT AS IN DOCUMENTATION & BELOW; IT'S NOT SUPPOSED TO BE "SIGN(GAMMAY)"

EX = (PEX1 + PEX2 .* DFZ +PEX3.*DFZ.*DFZ ) .* (1-PEX4 .*sign(KAPPAX)) .* LEX ; %24
% NOTE: LVX MULTIPLIES ONLY PVX1&2 IN DOCUMENTATION; ORIG VERSION MULT ALL TERMS

SVX = FZ .* (PVX1 + PVX2 .* DFZ) .* LVX .* LMUX .* ZETA1; %28

FX0 = DX .* sin((CX .* atan(BX .* KAPPAX -EX .*(BX.*KAPPAX - atan(BX.*KAPPAX))))+SVX);


%--Longitudinal Force combined slip---

SHXA = RHX1;

CXA = RCX1;

EXA = REX1 + REX2.* DFZ;

BXA = (RBX1 + RBX3.* GAMMAX.* GAMMAX).* cos(atan(RBX2.* KAPPA)).* LXA;

ALPHA_S = sin(ALPHA) + SHXA;

GXA0 = cos(CXA.*atan(BXA.*SHXA - EXA.*(BXA.*SHXA - atan(BXA.*SHXA))));

GXA = cos(CXA.*atan(BXA.*ALPHA_S - EXA.*(BXA.* ALPHA_S - atan(BXA.*ALPHA_S))))./ GXA0;

FX = GXA .* FX0;
