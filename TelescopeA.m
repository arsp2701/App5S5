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

%% Fonction de transfert de base
numAZ = [1.59e09];
denAZ = [1 1020.51 25082.705 3102480.725 64155612.5 8.27e7 0];
TFaz  = tf(numAZ,denAZ);
figure(1)
margin(TFaz)
figure(2)
rlocus(TFaz)

%% Avance de phase / PD en az ta
temps = linspace(0,20,1000)';
thetai = atan(-pi()/log(MP_A_az/100));
zetaaz = cos(thetai);
% zetaaz = 0.608334677523174;
% zetaaz = (zetaa+119*0.01);
wn1az = (4/(Ts_A_az*zetaaz));
wn2az = ((1+1.1*zetaaz+1.4*zetaaz^2)/Tm90_A_az);
wn3az = ((pi()-acos(zetaaz))/(Tm100_A_az*sqrt(1-zetaaz^2)));
% trouver le wn maximale
wnaz  = max([wn1az,wn2az,wn3az]);

% wnaz = wna+((35-1)*0.01);
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
c = 11;
a = 17;
alpha  = 180-(phivaz);
phizaz = (alpha+(dphiaz+a))/2;
phipaz =(alpha-(dphiaz+a))/2;
%Poles et zero voulu
Zaz = Praz-(Piaz/tand(phizaz));
Paz = Praz-(Piaz/tand(phipaz));
numZaz = [1 -Zaz];
denPaz = [1 -Paz];
TfZPaz = tf(numZaz,denPaz);

Kaaz = 1+(0.1*c)/norm(((polyval(numZaz,Poleaz)/polyval(denPaz,Poleaz))*(polyval(numAZ,Poleaz)/polyval(denAZ,Poleaz))));
[b,a] = cheby1(1,0.9,[40 75],'stop','s');
cp = tf(b,a);
% figure()
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(Kaaz*TfZPaz*TFaz)
figure()
margin(cp*Kaaz*TfZPaz*TFaz);
% figure()
% margin(Kaaz*TfZPaz*TFaz);
[GM,PM,Wp,Wg] = margin(cp*Kaaz*TfZPaz*TFaz);
PMaz = PMtaz*(180/pi)*Wg
PMtazr = (PM*(pi/180))/Wg
TFfazf = feedback(cp*Kaaz*TfZPaz*TFaz,1);
figure
lsim(TFfazf,ones(size(temps)),temps);
%step(TFfazf)
hold on 
plot(temps,1.02*ones(size(temps)),'-.b')
plot(temps,0.98*ones(size(temps)),'-.b')
plot(temps,1.3*ones(size(temps)),'-.b')
xline(1.25);
[numapaz,denapaz] = tfdata(cp*Kaaz*TfZPaz*TFaz,'v');
denapaz = denapaz(1:9);
Kvelaz = abs((polyval(numapaz,0)/polyval(denapaz,0)));
errper = 1/Kvelaz
figure()
test = lsim(TFfazf,temps,temps);
plot(temps,temps-test)
hold on
plot(temps,Err_r_A_az*1.02*ones(size(temps)),'-.b')
plot(temps,Err_r_A_az*0.98*ones(size(temps)),'-.b')
xline(1.25);

%trajectoire

Profile_Tracking
figure()
plot(ttrk,utrk)
title('Trajectoire réelle et suivie')
xlabel('Temps (s)')
ylabel('Position')
hold on 
trajectoire = lsim(TFfazf,utrk,ttrk);
plot(ttrk,trajectoire)
legend ('Trajectoire réelle','Trajectoire suivie')

figure ()
plot(ttrk,(utrk-trajectoire))
title('Erreur entre la trajectoire voulue et l''asservissement')
xlabel('Temps (s)')
ylabel('Erreur valeur absolue')

% figure ()
% plot(ttrk(1),(utrk-trajectoire)./utrk)
% ylim([-0.5 1])
% title('Erreur entre la trajectoire voulue et l''asservissement en %')
% xlabel('Temps (s)')
% ylabel('Erreur (%)')

