close all
clear all
clc
%  Téléscope A

%% Élément demandés antenne A
% Élément azimut
MP_A_az = 25;
Ts_A_az = 1;
Tm90_A_az = 0.15;
Tm100_A_az = 0.25;

% Réponse à l'échelon
Err_e_A_az = 0;
% Réponse à une rampe
Err_r_A_az = 0.05; %deg

% CRITÈRES DE SÉCURITÉ
MGtaz = 10; %dB
PMtaz = 0.1; %s
ATvraz = -15; %dB

%CRITÈRES FINAUX
MPtazac = 30;
% Réponse à l'échelon
Tstazeac = 1.25; %s
% Réponse à une rampe
Tstazrac = 1.25; %s

% Élément élévation
% Réponse à une rampe
Errrtae = 0;
% Réponse à une parabole
Errptae = 0.1; %deg

PMtae = 0.08; %s
ATvrae = -15; %dB
%Critères finaux
MPtazec = 35;
% Réponse à l'échelon
Tstazeec = 1.5; %s
% Réponse à une rampe
Tstazrec = 1.5; %s
% Réponse à une parabole
Tstazpec = 3; %s

%%
%% Fonction de transfert de base
numAZ = [1.59e09];
denAZ = [1 1020.51 25082.705 3102480.725 64155612.5 8.27e7 0];
TFaz  = tf(numAZ,denAZ);
figure(1)
margin(TFaz)
figure(2)
rlocus(TFaz)
numEL  = [7.95e09];
denEL  = [1 1020.51 37082.705 15346520.725 320776412.5 4.135e8 0];
TFe   = tf(numEL,denEL);
figure(3)
margin(TFe)
figure(4)
rlocus(TFe)
%% Avance de phase / PD en az ta
thetai = atan(-pi()/log(MP_A_az/100));
zetaaz = cos(thetai);
wn1az = 4/(Ts_A_az*zetaaz);
wn2az = (1+1.1*zetaaz+1.4*zetaaz^2)/Tm90_A_az;
wn3az = (pi()-acos(zetaaz))/(Tm100_A_az*sqrt(1-zetaaz^2));
% trouver le wn maximale
wnaz  = max([wn1az,wn2az,wn3az]);
Praz  = -zetaaz*wnaz;
Piaz  = wnaz*sqrt(1-zetaaz^2);
Poleaz= Praz+Piaz*1i;
% figure(5)
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(TFe)
%Angle des poles et zero voulu
dphiaz = 360+(-180-rad2deg(angle(polyval(numAZ,Poleaz)))+rad2deg(angle(polyval(denAZ,Poleaz))));
phivaz = 180-rad2deg(angle(Poleaz));
alpha  = 180-phivaz;
phizaz = (alpha+dphi)/2;
phipaz =(alpha-dphi)/2;
%Poles et zero voulu
Zaz = Praz-Piaz/tand(phizaz);
Paz = Praz-Piaz/tand(phipaz);
numZaz = [1 -Zaz];
denPaz = [1 -Paz];
TfZPaz = tf(numZaz,denPaz);
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(TfZPaz*TFaz)
% figure()
% margin(TfZPaz*TFaz)

