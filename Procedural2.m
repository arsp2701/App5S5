clear all
close all 
clc

% %%
% 
% num = [5 50];
% den = [1 2 2];
% Ts = 0.2;
% 
% 
% FTBO = tf(num,den);
% 
% [num_cap,den_cap] = pade(Ts,2);
% 
% FTBO_cap = tf(num_cap,den_cap);
% 
% rlocus(FTBO*FTBO_cap)
%% E2 procedural 2

num = [1 4];
den = [1 8 3 0];
FTBO = tf(num,den);
margin(FTBO)

PM = 50;
BW = 5.2;
err_ramp = 0.005;

zeta_d = 1/2*sqrt(tand(PM)*sind(PM));

Wg = BW*sqrt(sqrt(1+4*zeta_d^4)-2*zeta_d^2)/((1-2*zeta_d^2)+sqrt(4*zeta_d^4-4*zeta_d^2+2))^(1/2);
K_star = 1/abs(polyval(num,Wg*1i)/polyval(den,Wg*1i));

[GM_o,PM_o,Wp_o,Wg_o] = margin(K_star*FTBO);


delta_phi = PM - PM_o + 6;

alpha = (1-sind(delta_phi))/(1+sind(delta_phi));

T = 1/Wg/sqrt(alpha);

zero_c = -1/T;
pole_c = -1/(alpha*T);



Ka = K_star/sqrt(alpha);

num_c = [1 -zero_c];
den_c = [1 -pole_c];

Gc = tf(num_c,den_c);

FTBO_Av = Ka*Gc*FTBO;
% 
% bode(filtre)
% figure
% margin(FTBO_Av*filtre)

%% Retard de phase avec

% 0.005 err ramp


Kvel_v = 1/err_ramp;
Kvel_a = 4/3;

k_star = Kvel_v/Kvel_a;


beta = k_star;

zero_c_r = -Wg/10;
pole_c_r = -Wg/beta/10;

num_r = [1 -zero_c_r];
den_r = [1 -pole_c_r];

comp_r = tf(num_r,den_r);
Kr = k_star/(beta*abs((1i*Wg-zero_c_r)/(1i*Wg-pole_c_r)));

% margin(FTBO_Av*Kr*comp_r)




%% Méthode2

num = [1 4];
den = [1 8 3 0];

FTBO = tf(num,den);

err_v = 0.005;
PM_v = 50;

K_vel_v = 1/err_v;

K_vel = 4/3;

K_v = K_vel_v/K_vel;

% margin(FTBO*K_v)







