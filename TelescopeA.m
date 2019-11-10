close all
clear all
clc
%  Téléscope A

%% Élément demandés antenne A
% Élément azimut
MP_A_az = 25;
Ts_A_az = 1.25;
Tm90_A_az = 0.25;
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
close all
temps = linspace(0,2,1000)';
Ts_A_az = 0.107;

MP_A_az = 9;
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
phizaz = (alpha+dphiaz)/2;
phipaz =(alpha-dphiaz)/2;
%Poles et zero voulu
Zaz = Praz-(Piaz/tand(phizaz));
Paz = Praz-(Piaz/tand(phipaz));
numZaz = [1 -Zaz];
denPaz = [1 -Paz];
TfZPaz = tf(numZaz,denPaz);
Kaaz = 1/abs(((polyval(numZaz,Poleaz)/polyval(denPaz,Poleaz))*(polyval(numAZ,Poleaz)/polyval(denAZ,Poleaz))));
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(Kaaz*TfZPaz*TFaz)
% figure(6)
% margin(Kaaz*TfZPaz*TFaz);
%test = lsim(Kaaz*TfZPaz*TFaz,temps,temps);
TFfazf = feedback(Kaaz*TfZPaz*TFaz,1);
% figure(7)
% test = lsim(TFfazf,temps,temps);
% plot(temps,temps-test)
% figure()
% lsim(TFfazf,ones(size(temps)),temps);
% hold on 
% plot(temps,1.02*ones(size(temps)),'-.b')
% plot(temps,0.98*ones(size(temps)),'-.b')
% plot(temps,1.3*ones(size(temps)),'-.b')
[numapaz,denapaz] = tfdata(Kaaz*TfZPaz*TFaz,'v');
denapaz = denapaz(1:7);
Kvelaz = abs((polyval(numapaz,0)/polyval(denapaz,0)));
errper = 1/Kvelaz;
[GM,PM,Wp,Wg] = margin(Kaaz*TfZPaz*TFaz);
PMaz = PMtaz*(180/pi)*Wg;
% PI
Ketaz = (1/Err_r_A_az)/Kvelaz;

PZpiaz = Praz/11;
%PPpiaz = 0;
PPpiaz = PZpiaz/(Ketaz*10.3);
numpiZaz = [1 -PZpiaz];
denpiPaz = [1 -PPpiaz];
TfZPiaz = tf(numpiZaz,denpiPaz);
Kraz = 1/abs(((polyval(numpiZaz,Poleaz)/polyval(denpiPaz,Poleaz))*(polyval(numapaz,Poleaz)/polyval(denapaz,Poleaz))));
TFfazpi = feedback(Kaaz*TfZPaz*TFaz*Kraz*TfZPiaz,1);
% figure()
margin(Kaaz*TfZPaz*TFaz*Kraz*TfZPiaz) %pour le moment seul le pic pause probleme
% figure()
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(Kaaz*TfZPaz*TFaz*Kraz*TfZPiaz)
% lsim(TFfazpi,ones(size(temps)),temps);
% axis([1 2 0.95 1.03])

