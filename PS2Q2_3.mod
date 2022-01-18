// Neoclassical Growth Model
// Transitory L shock
// dynare tran_L_shock noclearall;

var K Y C I r w;
varexo e;
parameters cBETA cSIGMA cALPHA cDELTA ASS  LSS eSS eSS0;

cBETA = 0.995; 
cSIGMA = 2.0; 
cALPHA = 0.3; 
cDELTA= 0.25; 
ASS=1.0; 
LSS=1.0;
eSS=1.0; 
eSS0=0.5;

KSS = eSS*LSS*((cBETA^(-1)+cDELTA-1)/(cALPHA*ASS))^(1/(cALPHA-1));
CSS = ASS*KSS^cALPHA*(eSS*LSS)^(1-cALPHA) - cDELTA*KSS;
YSS = ASS*KSS^cALPHA*(eSS*LSS)^(1-cALPHA); 
rSS = ASS*cALPHA*KSS^(cALPHA-1)*(eSS*LSS)^(1-cALPHA); 
wSS = ASS*(1-cALPHA)*KSS^cALPHA*(eSS*LSS)^(-cALPHA); 
ISS = YSS-CSS;

model;
0 = C(+1)^cSIGMA/C^cSIGMA - cBETA*(ASS*cALPHA*K^(cALPHA-1)*(e*LSS)^(1-cALPHA) - cDELTA + 1);
0 = K - (ASS*K(-1)^cALPHA*(e*LSS)^(1-cALPHA) + (1-cDELTA)*K(-1) - C);
Y = ASS*K(-1)^cALPHA*(e*LSS)^(1-cALPHA); 
r = ASS*cALPHA*K(-1)^(cALPHA-1)*(e*LSS)^(1-cALPHA); 
w = ASS*(1-cALPHA)*K(-1)^cALPHA*(e*LSS)^(-cALPHA); 
I = Y-C;
end;

%================= Phase diagram
initval;
C = CSS;
K = 0.5*KSS;
Y = YSS;
r = rSS;
w = wSS;
I = ISS;
e = eSS;
end;

perfect_foresight_setup(periods=100);
perfect_foresight_solver;

IRFs_b1 = oo_.endo_simul;
IRFs_e1 = oo_.exo_simul;

initval;
C = CSS;
K = 2*KSS;
Y = YSS;
r = rSS;
w = wSS;
I = ISS;
e = eSS;
end;

perfect_foresight_setup(periods=100);
perfect_foresight_solver;

IRFs_b2 = oo_.endo_simul;
IRFs_e2 = oo_.exo_simul;

%================= IRFs

initval;
C = CSS;
K = KSS;
Y = YSS;
r = rSS;
w = wSS;
I = ISS;
e = eSS;
end;

shocks;
var e;
periods 1:5;
values 0.5;
end;

Tdur=5+1; // duration

perfect_foresight_setup(periods=100);
perfect_foresight_solver;

IRFs_b = oo_.endo_simul;
IRFs_e = oo_.exo_simul;

//var K Y C I r w;

Tfig=30;

figure(5);
subplot(3,3,1)
plot([0:Tfig],[KSS IRFs_b(1,1:Tfig)],'b','LineWidth',1.5);
title('K')
xlim([0 Tfig])
xline(Tdur)
subplot(3,3,2)
plot([0:Tfig],[YSS IRFs_b(2,2:Tfig+1)],'b','LineWidth',1.5);
title('Y')
xlim([0 Tfig])
xline(Tdur)
subplot(3,3,3)
plot([0:Tfig],[CSS IRFs_b(3,2:Tfig+1)],'b','LineWidth',1.5);
title('C')
xlim([0 Tfig])
xline(Tdur)
subplot(3,3,4)
plot([0:Tfig],[ISS IRFs_b(4,2:Tfig+1)],'b','LineWidth',1.5);
title('I')
xlim([0 Tfig])
xline(Tdur)
subplot(3,3,5)
plot([0:Tfig],[rSS IRFs_b(5,2:Tfig+1)],'b','LineWidth',1.5);
title('r')
xlim([0 Tfig])
xline(Tdur)
subplot(3,3,6)
plot([0:Tfig],[wSS IRFs_b(6,2:Tfig+1)],'b','LineWidth',1.5);
title('w')
xlim([0 Tfig])
xline(Tdur)
subplot(3,3,7)
plot([0:Tfig],[eSS IRFs_e(2:Tfig+1,1)'],'b','LineWidth',1.5);
title('L')
xlim([0 Tfig])
xline(Tdur)

Kval = 0.2:0.01:2;

figure(7);
plot([KSS IRFs_b(1,1:Tfig)],IRFs_b(3,1:Tfig+1),'r','LineWidth',2.5);
xlabel('K_t')
ylabel('C_t')
xlim([min(Kval) max(Kval)])
hold on
plot(IRFs_b1(1,1:Tfig),IRFs_b1(3,2:Tfig+1),'r--','LineWidth',1.5);
xlim([min(Kval) max(Kval)])
hold on
plot(IRFs_b2(1,1:Tfig),IRFs_b2(3,2:Tfig+1),'r--','LineWidth',1.5);
hold on
plot(Kval,ASS*Kval.^cALPHA*(eSS0*LSS)^(1-cALPHA)-cDELTA*Kval,'k--','LineWidth',1.5);
hold on
plot(Kval,ASS*Kval.^cALPHA*(eSS*LSS)^(1-cALPHA)-cDELTA*Kval,'k','LineWidth',1.5);
hold on
xline(eSS0*LSS*(1/ASS*1/cALPHA*(1/cBETA+cDELTA-1))^(1/(cALPHA-1)),'k--','LineWidth',1.5);
hold on
xline(eSS*LSS*(1/ASS*1/cALPHA*(1/cBETA+cDELTA-1))^(1/(cALPHA-1)),'k','LineWidth',1.5);
hold on

