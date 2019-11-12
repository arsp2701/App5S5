clear all
close all
clc

% Élément élévation
MP_A_ae = 25;
Ts_A_ae = 1;
Tm90_A_ae = 0.15;
Tm100_A_ae = 0.25;
% Réponse à une rampe
Errrtae = 0;
% Réponse à une parabole
Errptae = 0.1; %deg
MGtae = 10; %dB
PMtae = 0.08; %s
ATvrae = -15; %dB
%Critères finaux
MPtaec = 35;
% Réponse à l'échelon
Tstaeec = 1.5; %s
% Réponse à une rampe
Tstarec = 1.5; %s
% Réponse à une parabole
Tstapec = 3; %s

numae  = [7.95e09];
denae  = [1 1020.51 37082.705 15346520.725 320776412.5 4.135e8 0];
TFae   = tf(numae,denae);
figure(3)
margin(TFae)
figure(4)
rlocus(TFae)

%% Avance de phase / PD en ae 
temps = linspace(0,10,100000)';
thetai = atan(-pi()/log(MP_A_ae/100));
zetaae = cos(thetai);
wn1ae = (4/(Ts_A_ae*zetaae));
wn2ae = ((1+1.1*zetaae+1.4*zetaae^2)/Tm90_A_ae);
wn3ae = ((pi()-acos(zetaae))/(Tm100_A_ae*sqrt(1-zetaae^2)));
% trouver le wn maximale
wnae  = max([wn1ae,wn2ae,wn3ae]);
Prae  = -zetaae*wnae;
Piae  = wnae*sqrt(1-zetaae^2);
Poleae= Prae+Piae*1i;
% figure(5)
% plot(Prae,Piae,'p')
% hold on
% plot(Prae,-Piae,'p')
% rlocus(TFae)
%Angle des poles et zero voulu
dphiae = 360+(-180-rad2deg(angle(polyval(numae,Poleae)))+rad2deg(angle(polyval(denae,Poleae))));
phivae = 180-rad2deg(angle(Poleae));
alpha  = 180-(phivae);

phizae = (alpha+(dphiae+17))/2; % 13
phipae =(alpha-(dphiae+17))/2;
%Poles et zero voulu
Zae = Prae-(Piae/tand(phizae));
Pae = Prae-(Piae/tand(phipae));
numZae = [1 -Zae];
denPae = [1 -Pae];
TfZPae = tf(numZae,denPae);
Kaae = (1.1+(5.35*0.01))/norm(((polyval(numZae,Poleae)/polyval(denPae,Poleae))*(polyval(numae,Poleae)/polyval(denae,Poleae))));
TFaf = (Kaae*TfZPae*TFae);
[numiae deniae] = tfdata(TFaf,'v');
% figure()
% plot(Prae,Piae,'p')
% hold on
% plot(Prae,-Piae,'p')
% rlocus(Kaae*TfZPae*TFae)
% figure()
% margin(Kaae*TfZPae*TFae);
[GM,PM,Wp,Wg] = margin(Kaae*TfZPae*TFae);
PMae = PMtae*(180/pi)*Wg;
TFfaef = feedback(Kaae*TfZPae*TFae,1);
% figure()
% lsim(TFfaef,ones(size(temps)),temps);
% %step(TFfaef)
% hold on 
% plot(temps,1.02*ones(size(temps)),'-.b')
% plot(temps,0.98*ones(size(temps)),'-.b')
% plot(temps,1.35*ones(size(temps)),'-.b')
% xline(1.25);
% figure()
% test = lsim(TFfaef,temps,temps);
% plot(temps,temps-test)
% hold on
% plot(temps,Errptae*1.02*ones(size(temps)),'-.b')
% plot(temps,Errptae*0.98*ones(size(temps)),'-.b')
% xline(1.25);
% PI


PZpiae = Prae/4;
PPpiae = 0;
numpiZae = [1 -PZpiae];
denpiPae = [1 -PPpiae];
TfZPiae = tf(numpiZae,denpiPae);

Krae = 1/norm(((polyval(numpiZae,Poleae)/polyval(denpiPae,Poleae))*(polyval(numiae,Poleae)/polyval(deniae,Poleae))));
[b,a] = cheby1(1,1,[115 130],'stop','s');
cp = tf(b,a);
TFfaepi = feedback(cp*Kaae*TfZPae*TFae*Krae*TfZPiae,1);
figure()
margin(cp*Kaae*TfZPae*TFae*Krae*TfZPiae) 
figure()
plot(Prae,Piae,'p')
hold on
plot(Prae,-Piae,'p')
rlocus(cp*Kaae*TfZPae*TFae*Krae*TfZPiae)
figure()
test = lsim(TFfaepi,ones(size(temps)),temps);
%step(TFfaepi)
plot(temps,test)
hold on 
plot(temps,1.02*ones(size(temps)),'-.b')
plot(temps,0.98*ones(size(temps)),'-.b')
plot(temps,1.35*ones(size(temps)),'-.b')
xline(1.25);
[GM,PM,Wp,Wg] = margin(Kaae*TfZPae*TFae*Krae*TfZPiae);
PMae = PMtae*(180/pi)*Wg
PMtaer = PM*(pi/180)/Wg
figure()
test1 = lsim(TFfaepi,temps,temps);
plot(temps,temps-test1)
[numapae,denapae] = tfdata(Kaae*TfZPae*TFae*Krae*TfZPiae,'v');
denapae = denapae(1:7);
Kaccae = abs((polyval(numapae,0)/polyval(denapae,0)));
errper = 1/Kaccae
figure()
test2 = lsim(TFfaepi,0.5*temps.^2,temps);
plot(temps,((0.5*temps.^2)-test2));
hold on
plot(temps,1.02*errper*ones(size(temps)),'-.b')
plot(temps,0.98*errper*ones(size(temps)),'-.b')
xline(3);

%trajectoire

Profile_Tracking
figure(8)
plot(ttrk,utrk)
title('Trajectoire réelle et suivie')
xlabel('Temps (s)')
ylabel('Position')
hold on 
trajectoire = lsim(TFfaepi,utrk,ttrk);
plot(ttrk,trajectoire)
legend ('Trajectoire réelle','Trajectoire suivie')

figure (9)
plot(ttrk,(utrk-trajectoire))
title('Erreur entre la trajectoire voulue et l''asservissement')
xlabel('Temps (s)')
ylabel('Erreur valeur absolue')

figure (10)
plot(ttrk,(utrk-trajectoire)./utrk)
ylim([-0.5 1])
title('Erreur entre la trajectoire voulue et l''asservissement en %')
xlabel('Temps (s)')
ylabel('Erreur (%)')