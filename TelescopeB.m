close all
clear all
clc
% Téléscope B

%% Éléments demandés antenne B
BW = 10; %rad/s
PM = 50; % deg +-5
% Réponse à une rampe
err_r = 0.005;


% CRITÈRES DE SÉCURITÉ
GMbz = 15; %dB
ATTvbz = -12; %dB
% B7 à calculer
%Critère finaux
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

delta_phi = PM - PM_o+18;

alpha = (1-sind(delta_phi))/(1+sind(delta_phi));

T = 1/Wg/sqrt(alpha);

zero_c = -1/T;
pole_c = -1/(alpha*T);

Ka = K_star/sqrt(alpha);

num_c = [1 -zero_c];
den_c = [1 -pole_c];

Gc = tf(num_c,den_c);

FTBO_Av = Ka*Gc*FTBO;
% figure
% margin(FTBO_Av)

% retard de phase
Kvel_star = 1/err_r;

[num_av,den_av] = tfdata(FTBO_Av,'v');

den_v = den_av(1:end-1);

K_star_r = Kvel_star/(polyval(num_av,0)/polyval(den_v,0));

beta = K_star_r;

zero_re = -Wg/10;
pole_re = -Wg/10/beta;

retard = tf([1 -zero_re],[1 -pole_re]);

FTBO_Av_Re = FTBO_Av*retard;
% 
% figure(69)
%  margin(FTBO_Av_Re)

[num_F, den_F] = butter(1,[50 60],'stop','s');

filter = tf(num_F,den_F);
figure(1)
margin(FTBO_Av_Re*filter)
% 
x = linspace(0,15,10000)';
% x2 = 1.02*x;
% x3 = .98*x;
% % 
feed = feedback(FTBO_Av_Re*filter,1);
y = lsim(feed,x,x);
z = x-y;
figure(2)
plot(x,z);
title('Erreur en régime permanent à une rampe unitaire, B, Azimut')
% hold on
% % plot(x,x2)
% % plot(x,x3)
% legend

% calcul de la marge restante 
w5 = 1.02962; % vu sur le margin
[mag,phase] = bode(FTBO_Av_Re*filter,w5);
marge_phase = phase-180;
retard_max = marge_phase*pi/180/w5;

[a,b,c,d]=margin(FTBO_Av_Re*filter);

%% Fonction de transfert de base en élévation
numEL  = [7.95e09];
denEL  = [1 1020.51 37082.705 15346520.725 320776412.5 4.135e8 0];
TFe   = tf(numEL,denEL);
% figure(3)
% margin(TFe)

% figure(1)
% margin(TFaz)
% figure(2)
% rlocus(TFaz)

FTBO = tf(numEL,denEL);

zeta_d = 1/2*sqrt(tand(PM)*sind(PM));

Wg = BW*sqrt(sqrt(1+4*zeta_d^4)-2*zeta_d^2)/((1-2*zeta_d^2)+sqrt(4*zeta_d^4-4*zeta_d^2+2))^(1/2);
K_star = 1/abs(polyval(numEL,Wg*1i)/polyval(denEL,Wg*1i));

[GM_o,PM_o,Wp_o,Wg_o] = margin(K_star*FTBO);

delta_phi = PM - PM_o+10;

alpha = (1-sind(delta_phi))/(1+sind(delta_phi));

T = 1/Wg/sqrt(alpha);

zero_c = -1/T;
pole_c = -1/(alpha*T);

Ka = K_star/sqrt(alpha);

num_c = [1 -zero_c];
den_c = [1 -pole_c];

Gc = tf(num_c,den_c);

FTBO_Av = Ka*Gc*FTBO;

% retard de phase
Kvel_star = 1/err_r;

[num_av,den_av] = tfdata(FTBO_Av,'v');

den_v = den_av(1:end-1);

K_star_r = Kvel_star/(polyval(num_av,0)/polyval(den_v,0));

beta = K_star_r;

zero_re = -Wg/10;
pole_re = -Wg/10/beta;

retard = tf([1 -zero_re],[1 -pole_re]);

FTBO_Av_Re = FTBO_Av*retard;


[num_F, den_F] = butter(1,[113 133],'stop','s');

filter = tf(num_F,den_F);




figure(3)
margin(FTBO_Av_Re*filter)

x = linspace(0,15,10000)';
% x2 = 1.02*x;
% x3 = .98*x;
% % 
feed = feedback(FTBO_Av_Re*filter,1);
y = lsim(feed,x,x);
z = x-y;
figure(4)
plot(x,z);
% hold on
% plot(x,x2)
% plot(x,x3)

% calcul de la marge restante 
w5 = 1.3075; % vu sur le margin et vérifié avec bode w5
[mag,phase] = bode(FTBO_Av_Re*filter,w5);
marge_phase = phase-180;
retard_max = marge_phase*pi/180/w5;



