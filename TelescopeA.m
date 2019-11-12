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
temps = linspace(0,2,1000)';
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
%  for c = 0.1:.1:0.2
% %    figure()
% disp(['--------------------',num2str(c),'--------------------'])
%      for ac = 12:.1:18
c = 0;
ac = 2;
alpha  = 180-(phivaz);
phizaz = (alpha+(dphiaz+ac))/2;
phipaz =(alpha-(dphiaz+ac))/2;
%Poles et zero voulu
Zaz = Praz-(Piaz/tand(phizaz));
Paz = Praz-(Piaz/tand(phipaz));
numZaz = [1 -Zaz];
denPaz = [1 -Paz];
TfZPaz = tf(numZaz,denPaz);

Kaaz = (1+(0.1*c))/norm(((polyval(numZaz,Poleaz)/polyval(denPaz,Poleaz))*(polyval(numAZ,Poleaz)/polyval(denAZ,Poleaz))));
% figure()
% plot(Praz,Piaz,'p')
% hold on
% plot(Praz,-Piaz,'p')
% rlocus(Kaaz*TfZPaz*TFaz)
% figure()
% margin(cp*Kaaz*TfZPaz*TFaz);
% figure()
% margin(Kaaz*TfZPaz*TFaz);
[GM,PM,Wp,Wg] = margin(Kaaz*TfZPaz*TFaz);
PMaz = PMtaz*(180/pi)*Wg;
PMtazr = (PM*(pi/180))/Wg;
TFfazf = feedback(Kaaz*TfZPaz*TFaz,1);
% figure()
test1 = lsim(TFfazf,ones(size(temps)),temps);
maxp = max(test1);
% %step(TFfazf)
% hold on 
% plot(temps,1.02*ones(size(temps)),'-.b')
% plot(temps,0.98*ones(size(temps)),'-.b')
% plot(temps,1.3*ones(size(temps)),'-.b')
% xline(1.25);
% title('step ap')
[numapaz,denapaz] = tfdata(Kaaz*TfZPaz*TFaz,'v');
denapaz = denapaz(1:7);
Kvelaz = abs((polyval(numapaz,0)/polyval(denapaz,0)));
errper = 1/Kvelaz;

% figure()
% test = lsim(TFfazf,temps,temps);
% plot(temps,temps-test)
% hold on
% plot(temps,Err_r_A_az*1.02*ones(size(temps)),'-.b')
% plot(temps,Err_r_A_az*0.98*ones(size(temps)),'-.b')
% xline(1.25);
% title('ramp ap')
disp(['----------',num2str(ac),'----------'])
disp(['GM = ',num2str(20*log10(GM))])
disp(['Pm des = ', num2str(PMaz)])
disp(['PM = ',num2str(PM)])
disp(['marge lousse = ',num2str(PM-PMaz)])
disp(['maxp = ',num2str(maxp)])
disp(['errper = ',num2str(errper)])
%      end 
%  end
%% avance 2

[b,a] = cheby1(2,1,[43 64.5],'stop','s');
cp = tf(b,a);
PMet = PMaz;
[numaf1 denaf1] =tfdata(Kaaz*TfZPaz*TFaz,'v');
Ket = 1/norm(polyval(numaf1,Wg*1i)/polyval(denaf1,Wg*1i));

Dphi2 = PMet-PM;
alpha2 = (1-sind(Dphi2))/(1+sind(Dphi2));
T = 1/(Wg*sqrt(alpha2));
Zaz2 = -1/(T);
Paz2 = -1/(alpha2*T);
numZaz2 = [1 -Zaz2];
denPaz2 = [1 -Paz2];
Ka2 = Ket/sqrt(alpha2);
TfZPaz2 = 0.98*Ka2*tf(numZaz2,denPaz2);
TFfaz2  = (cp*TfZPaz2*Kaaz*TfZPaz*TFaz);
TFfazf2 = feedback(TFfaz2,1);
% figure()
margin(TFfaz2)
hold on
[GM2,PM2,Wp2,Wg2] = margin(TFfaz2);
PMaz2 = PMtaz*(180/pi)*Wg2;
PMtazr2 = (PM2*(pi/180))/Wg2;

[numapaz2,denapaz2] = tfdata(TFfaz2,'v');
denapaz2 = denapaz2(1:12);
Kvelaz = abs((polyval(numapaz,0)/polyval(denapaz2,0)));
errper = 1/Kvelaz

% figure()
% lsim(TFfazf2,ones(size(temps)),temps);
% hold on 
% plot(temps,1.02*ones(size(temps)),'-.b')
% plot(temps,0.98*ones(size(temps)),'-.b')
% plot(temps,1.3*ones(size(temps)),'-.b')
% xline(1.25);
% title('step ap')
% figure()
% test = lsim(TFfazf2,temps,temps);
% plot(temps,temps-test)
% hold on
% plot(temps,Err_r_A_az*1.02*ones(size(temps)),'-.b')
% plot(temps,Err_r_A_az*0.98*ones(size(temps)),'-.b')
% xline(1.25);
% title('ramp ap')
%% trajectoire

Profile_Tracking
figure()
plot(ttrk,utrk)
title('Trajectoire réelle et suivie')
xlabel('Temps (s)')
ylabel('Position')
hold on 
trajectoire = lsim(TFfazf2,utrk,ttrk);
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
