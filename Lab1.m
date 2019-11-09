clear all
clc
close all

%% E1 - Lab 

Mp = 0.06;
tm_10_90 = 0.004;
tp = 0.008;
ts = 0.010;
err = 0.00005;

num = 4500;
den = [1 361.2 0];

FTBO = tf(num,den);

phi = atand(1/log(Mp)*(-pi));
zeta = cosd(phi);

Wn1 = 4/(zeta*ts);
Wn2 = (1+1.1*zeta+1.4*zeta^2)/tm_10_90;
Wn3 = pi/(tp*(sqrt(1-zeta^2)));
% le Wn le plus grand est Wn1


poles = [-zeta*Wn1+1i*Wn1*sqrt(1-zeta^2) -zeta*Wn1-1i*Wn1*sqrt(1-zeta^2)];


%  plot(real(poles),imag(poles),'p')
%  hold on
% % rlocus(FTBO)

delta_phi = rad2deg(-pi - angle(polyval(num,poles(1))) + angle(polyval(den,poles(1))));

alpha = 180-phi;

phi_z = (alpha + delta_phi)/2;
phi_p = (alpha - delta_phi)/2;

zero_v = real(poles(1))-imag(poles(1))/tand(phi_z);
pole_v = real(poles(1))-imag(poles(1))/tand(phi_p);

num_c = [1 -zero_v];
den_c = [1 -pole_v];

Gc = tf (num_c,den_c);

Gs_v = polyval(num,poles(1))/polyval(den,poles(1));
Gc_v = (poles(1)-zero_v)/(poles(1)-pole_v);


Ka = 1/abs(Gc_v*Gs_v);
% 

% rlocus(Gc*FTBO)


%% c)

GS = Ka*Gc*FTBO;


[num_gs,den_gs] = tfdata(GS,'v');
Kvel_d = 1/err;

% car lim s-»0 = s*G(s)
Kvel = num_gs(4)/den_gs(3);
facteur = Kvel_d/Kvel;


FTBF = feedback(GS,1);



zero_retard = real(poles(1))/10;
pole_retard = zero_retard/facteur;


num_retard = [1 -zero_retard];
den_retard = [1 -pole_retard];

G_retard = tf(num_retard,den_retard);
% figure(1)
% rlocus(G_retard*GS)


%% E1-e)
%PD

zero_pd = real(poles(1))-imag(poles(1))/tand(delta_phi);


Gc_pd = (poles(1)-zero_pd);

K_pd = 1/abs(Gs_v*Gc_pd);

GC_pd = tf([1 -zero_pd],[1]);


% rlocus(K_pd*GC_pd*FTBO);

%PI


zero_pi = real(poles(1))/10;
Gc_pi = (poles(1)-zero_pi)/poles(1);
kp_pi = 1/abs(Gs_v*Gc_pi);
PI = tf([1 -zero_pi],[1 0]);

% figure(2)
% plot(real(poles),imag(poles),'p')
%  hold on
% rlocus(kp_pi*PI*FTBO*GC_pd*K_pd)


%% E2


num = [1 1];
den = conv([1 6 10],[1 7]);

FTBO = tf(num,den);

% rlocus(FTBO)



num_i = num;
den_i = conv([1 6 10],[1 7 0]);

FTBO_i = tf(num_i,den_i);




%Conception d'un PD qui a une avance de phase de delta_phi sur 2

ts_2 = 0.1;
MP = 0.05;

phi = atand(1/log(MP)*(-pi));
zeta = cosd(phi);
Wn = 4/zeta/ts_2;

poles = [(-zeta*Wn+1i*Wn*sqrt(1-zeta^2)),(-zeta*Wn-1i*Wn*sqrt(1-zeta^2))];
 
plot(real(poles),imag(poles),'p')
hold on
% rlocus(FTBO_i);

phase_FTBO_i = rad2deg(angle(polyval(num_i,poles(1)))-(angle(polyval(den_i,poles(1)))));

delta_phi = -180 - phase_FTBO_i;


delta_phi_2 = delta_phi/2;

zero_pd = real(poles(1))+imag(poles(1))/tan(delta_phi_2);


PD_1 = tf([1 -zero_pd],1);
Gs = (polyval(num_i,poles(1)))/polyval(den_i,poles(1));
Kd = 1/abs((poles(1)-zero_pd)^2*Gs);

rlocus(PD_1^2*Kd*FTBO_i)



