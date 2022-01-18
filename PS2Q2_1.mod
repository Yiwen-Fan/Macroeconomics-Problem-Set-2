// Neoclassical Growth Model
// Permanent L shock
// dynare Per_L_shock noclearall;

var K Y C I r w;
varexo e;
parameters cBETA cSIGMA cALPHA cDELTA ASS  LSS eSS eSS0;

cBETA = 0.995; 
cSIGMA = 2.0; 
cALPHA = 0.3; 
cDELTA= 0.25; 
ASS=1.0; 
LSS=1.0;
eSS0=1.0; 
eSS=0.5;

KSS0 = eSS0*LSS*((cBETA^(-1)+cDELTA-1)/(cALPHA*ASS))^(1/(cALPHA-1));
CSS0 = ASS*KSS0^cALPHA*(eSS0*LSS)^(1-cALPHA) - cDELTA*KSS0;
YSS0 = ASS*KSS0^cALPHA*(eSS0*LSS)^(1-cALPHA); 
rSS0 = ASS*cALPHA*KSS0^(cALPHA-1)*(eSS0*LSS)^(1-cALPHA); 
wSS0 = ASS*(1-cALPHA)*KSS0^cALPHA*(eSS0*LSS)^(-cALPHA); 
ISS0 = YSS0-CSS0;

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



initval;
C = CSS0;
K = KSS0;
Y = YSS0;
r = rSS0;
w = wSS0;
I = ISS0;
e = eSS0;
end;

endval;
C = CSS;
K = KSS;
Y = YSS;
r = rSS;
w = wSS;
I = ISS;
e = eSS;
end;

perfect_foresight_setup(periods=100);
perfect_foresight_solver;

IRFs_b = oo_.endo_simul;
IRFs_e = oo_.exo_simul;

//var K Y C I r w;

Tfig=30;

figure(5);
subplot(3,3,1)
plot([0:Tfig],[KSS0 IRFs_b(1,1:Tfig)],'b','LineWidth',1.5);
title('K')
xlim([0 Tfig])
subplot(3,3,2)
plot([0:Tfig],[YSS0 IRFs_b(2,2:Tfig+1)],'b','LineWidth',1.5);
title('Y')
xlim([0 Tfig])
subplot(3,3,3)
plot([0:Tfig],[CSS0 IRFs_b(3,2:Tfig+1)],'b','LineWidth',1.5);
title('C')
xlim([0 Tfig])
subplot(3,3,4)
plot([0:Tfig],[ISS0 IRFs_b(4,2:Tfig+1)],'b','LineWidth',1.5);
title('I')
xlim([0 Tfig])
subplot(3,3,5)
plot([0:Tfig],[rSS0 IRFs_b(5,2:Tfig+1)],'b','LineWidth',1.5);
title('r')
xlim([0 Tfig])
subplot(3,3,6)
plot([0:Tfig],[wSS0 IRFs_b(6,2:Tfig+1)],'b','LineWidth',1.5);
title('w')
xlim([0 Tfig])
subplot(3,3,7)
plot([0:Tfig],[eSS0 IRFs_e(2:Tfig+1,1)'],'b','LineWidth',1.5);
title('L')
xlim([0 Tfig])

Kval = 0.2:0.01:2;

figure(6);
plot([KSS0 IRFs_b(1,1:Tfig)],IRFs_b(3,1:Tfig+1),'r','LineWidth',2.5);
xlabel('K_t')
ylabel('C_t')
xlim([min(Kval) max(Kval)])
hold on
plot(IRFs_b1(1,1:Tfig),IRFs_b1(3,2:Tfig+1),'r--','LineWidth',1.5);
xlim([min(Kval) max(Kval)])
hold on
plot(IRFs_b2(1,1:Tfig),IRFs_b2(3,2:Tfig+1),'r--','LineWidth',1.5);
hold on
plot(Kval,ASS*Kval.^cALPHA*(eSS*LSS)^(1-cALPHA)-cDELTA*Kval,'k--','LineWidth',1.5);
hold on
plot(Kval,ASS*Kval.^cALPHA*(eSS0*LSS)^(1-cALPHA)-cDELTA*Kval,'k','LineWidth',1.5);
hold on
xline(eSS*LSS*(1/ASS*1/cALPHA*(1/cBETA+cDELTA-1))^(1/(cALPHA-1)),'k--','LineWidth',1.5);
hold on
xline(eSS0*LSS*(1/ASS*1/cALPHA*(1/cBETA+cDELTA-1))^(1/(cALPHA-1)),'k','LineWidth',1.5);
hold on

