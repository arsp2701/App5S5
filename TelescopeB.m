close all
clear all
clc
% T�l�scope B

%% �l�ments demand�s antenne B
BW = 10; %rad/s
PM = 50; % deg +-5
% R�ponse � une rampe
err_r = 0.005;


% CRIT�RES DE S�CURIT�
GMbz = 15; %dB
ATTvbz = -12; %dB
% B7 � calculer
%Crit�re finaux
Tsbz = 14; %s

%% Fonction de transfert de base azimut
numAZ = [1.59e09];
denAZ = [1 1020.51 25082.705 3102480.725 64155612.5 8.27e7 0];
TFaz  = tf(numAZ,denAZ);
% figure(1)
% margin(TFaz)
% figure(2)
% rlocus(TFaz)

FTBO = tf(numAZ,denAZ);

zeta_d = 1/2*sqrt(tand(PM)*sind(PM));

Wg = BW*sqrt(sqrt(1+4*zeta_d^4)-2*zeta_d^2)/((1-2*zeta_d^2)+sqrt(4*zeta_d^4-4*zeta_d^2+2))^(1/2);
K_star = 1/abs(polyval(numAZ,Wg*1i)/polyval(denAZ,Wg*1i));

[GM_o,PM_o,Wp_o,Wg_o] = margin(K_star*FTBO);

delta_phi = PM - PM_o+5;

alpha = (1-sind(delta_phi))/(1+sind(delta_phi));

T = 1/Wg/sqrt(alpha);

zero_c = -1/T;
pole_c = -1/(alpha*T);

Ka = K_star/sqrt(alpha)*.95;

num_c = [1 -zero_c];
den_c = [1 -pole_c];

Gc = tf(num_c,den_c);

FTBO_Av = Ka*Gc*FTBO;
% figure
% margin(FTBO_Av)
% figure
% rlocus(FTBO_Av)

[Wn,zeta] = damp(FTBO_Av)


%% Fonction de transfert de base en �l�vation
numEL  = [7.95e09];
denEL  = [1 1020.51 37082.705 15346520.725 320776412.5 4.135e8 0];
TFe   = tf(numEL,denEL);
% figure(3)
% margin(TFe)
% figure(4)
% rlocus(TFe)