close all
clear all
clc
%  T�l�scope A

%% �l�ment demand�s antenne A
% �l�ment azimut
MP_A_az = 25;
Ts_A_az = 1;
Tm90_A_az = 0.15;
Tm100_A_az = 0.25;

% R�ponse � l'�chelon
Err_e_A_az = 0;
% R�ponse � une rampe
Err_r_A_az = 0.05; %deg

% CRIT�RES DE S�CURIT�
MGtaz = 10; %dB
PMtaz = 0.1; %s
ATvraz = -15; %dB
%CRIT�RES FINAUX
MPtazac = 30;
% R�ponse � l'�chelon
Tstazeac = 1.25; %s
% R�ponse � une rampe
Tstazrac = 1.25; %s

% �l�ment �l�vation
% R�ponse � une rampe
Errrtae = 0;
% R�ponse � une parabole
Errptae = 0.1; %deg

PMtae = 0.08; %s
ATvrae = -15; %dB
%Crit�res finaux
MPtazec = 35;
% R�ponse � l'�chelon
Tstazeec = 1.5; %s
% R�ponse � une rampe
Tstazrec = 1.5; %s
% R�ponse � une parabole
Tstazpec = 3; %s

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
figure(5)
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(TFe)
%Angle des poles et zero voulu
dphiaz = 360+(-180-rad2deg(angle(polyval(numAZ,Poleaz)))+rad2deg(angle(polyval(denAZ,Poleaz))));
phivaz = 180-rad2deg(angle(Poleaz));
alpha  = 180-phivaz;
phizaz = (alpha+dphiaz)/2;
phipaz =(alpha-dphiaz)/2;
%Poles et zero voulu
Zaz = Praz-(Piaz/tand(phizaz));
Paz = Praz-(Piaz/tand(phipaz));
numZaz = [1 -Zaz];
denPaz = [1 -Paz];
TfZPaz = tf(numZaz,denPaz);
Kaaz = 1/abs(((polyval(numZaz,Poleaz)/polyval(denPaz,Poleaz))*(polyval(numAZ,Poleaz)/polyval(denAZ,Poleaz))));
plot(Praz,Piaz,'p')
hold on
plot(Praz,-Piaz,'p')
rlocus(Kaaz*TfZPaz*TFaz)
figure(6)
margin(Kaaz*TfZPaz*TFaz);
temps = linspace(0,2,100000);
%test = lsim(Kaaz*TfZPaz*TFaz,temps,temps);
TFfazf = feedback(Kaaz*TfZPaz*TFaz,1);
figure(7)
lsim(TFfazf,ones(size(temps)),temps);
hold on 
plot(temps,1.02*ones(size(temps)),'-.b')
plot(temps,0.98*ones(size(temps)),'-.b')
plot(temps,1.25*ones(size(temps)),'-.b')
[numapaz,denapaz] = tfdata(Kaaz*TfZPaz*TFaz,'v');
denapaz = denapaz(1:7);
Kvelaz = abs((polyval(numapaz,0)/polyval(denapaz,0)));
errper = 1/Kvelaz;
[GM,PM,Wp,Wg] = margin(Kaaz*TfZPaz*TFaz);
PMaz = PMtaz*(180/pi)*Wg;
%% PI
Ketaz = (1/Err_r_A_az)/Kvelaz;
PZpiaz = Praz/200;
PPpiaz = PZpiaz/Ketaz;
numpiZaz = [1 -PZpiaz];
denpiPaz = [1 -PPpiaz];
TfZPiaz = tf(numpiZaz,denpiPaz);
Kraz = 1/abs(((polyval(numpiZaz,Poleaz)/polyval(denpiPaz,Poleaz))*(polyval(numapaz,Poleaz)/polyval(denapaz,Poleaz))));
TFfazpi = feedback(Kaaz*TfZPaz*TFaz*Kraz*TfZPiaz,1);
figure()
margin(Kaaz*TfZPaz*TFaz*Kraz*TfZPiaz)
figure()
plot(Praz,Piaz,'p')
hold on
plot(Praz,-Piaz,'p')
rlocus(Kaaz*TfZPaz*TFaz*Kraz*TfZPiaz)
figure()
step(TFfazpi)