close all;
T = 50;
time=0:T;
cALPHA = 0.4;
cDELTA = 0.1;
cS=0.2;

K   = zeros(1,T);
Y   = zeros(1,T);
r   = zeros(1,T);
w   = zeros(1,T);
I   = zeros(1,T);
C   = zeros(1,T);


L   = ones(1,T);
A   = ones(1,T);


% A transitory drop in L from time 1 to 5
L(1:5) = 0.5;

Lss = 1;
Ass = 1;

% Steady states
Kss = (cS*Ass/cDELTA)^(1/(1-cALPHA))*Lss;
Yss = Ass*Kss^cALPHA*Lss^(1-cALPHA); 
rss = Ass*cALPHA*Kss^(cALPHA-1)*Lss^(1-cALPHA); 
wss = Ass*(1-cALPHA)*Kss^cALPHA*Lss^(-cALPHA); 
Iss = cS*Yss;
Css=Yss-Iss;


% Initial Condition
t=1;
K(t) = Kss;
Y(t) = A(t)*K(t).^cALPHA*L(t).^(1-cALPHA); 
r(t) = A(t)*cALPHA*K(t).^(cALPHA-1)*L(t).^(1-cALPHA); 
w(t) = A(t)*(1-cALPHA)*K(t).^cALPHA*L(t).^(-cALPHA); 
I(t) = cS*Y(t);
C(t) = Y(t)-I(t);
        
for t=2:T
	K(t) = I(t-1) + (1-cDELTA)*K(t-1);
	Y(t) = A(t)*K(t).^cALPHA*L(t).^(1-cALPHA); 
	r(t) = A(t)*cALPHA*K(t)^(cALPHA-1)*L(t).^(1-cALPHA); 
	w(t) = A(t)*(1-cALPHA)*K(t).^cALPHA*L(t).^(-cALPHA); 
	I(t) = cS*Y(t);
    C(t) = Y(t)-I(t);
end

subplot(2,4,1)
plot(time,[Kss K],'b','LineWidth',1.5);
hold on
plot(time,[Kss Kss*ones(1,T)],'r--');
hold on
xlim([0 T])
xlabel('Time')
title('Capital','FontSize',14)

subplot(2,4,2)
plot(time,[Iss I],'b','LineWidth',1.5);
hold on
plot(time,[Iss Iss*ones(1,T)],'r--');
hold on
xlim([0 T])
title('Investment','FontSize',14)
xlabel('Time')

subplot(2,4,3)
plot(time,[Css C],'b','LineWidth',1.5);
hold on
plot(time,[Css Css*ones(1,T)],'r--');
hold on
xlim([0 T])
title('Consumption','FontSize',14)
xlabel('Time')

subplot(2,4,4)
plot(time,[Yss Y],'b','LineWidth',1.5);
hold on
plot(time,[Yss Yss*ones(1,T)],'r--');
hold on
xlim([0 T])
xlabel('Time')
title('Output','FontSize',14)

subplot(2,4,5)
plot(time,[rss r],'b','LineWidth',1.5);
hold on
plot(time,[rss rss*ones(1,T)],'r--');
hold on
title('Rental price','FontSize',14)
xlim([0 T])
xlabel('Time')
subplot(2,4,6)
plot(time,[wss w],'b','LineWidth',1.5);
hold on
plot(time,[wss wss*ones(1,T)],'r--');
hold on
title('Wage','FontSize',14)
xlabel('Time')
xlim([0 T])
subplot(2,4,7)
plot(time,[Lss L],'b','LineWidth',1.5);
hold on
plot(time,[Lss Lss*ones(1,T)],'r--');
hold on
title('L','FontSize',14)
xlabel('Time')