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

%% Avance de phase / PD en ae ta
temps = linspace(0,2,1000)';
thetai = atan(-pi()/log(MP_A_ae/100));
zetaae = cos(thetai);
% zetaae = 0.608334677523174;
% zetaae = (zetaa+119*0.01);
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
% rlocus(TFe)
%Angle des poles et zero voulu
dphiae = 360+(-180-rad2deg(angle(polyval(numae,Poleae)))+rad2deg(angle(polyval(denae,Poleae))));
phivae = 180-rad2deg(angle(Poleae));
alpha  = 180-(phivae);
phizae = (alpha+(dphiae))/2;
phipae =(alpha-(dphiae))/2;
%Poles et zero voulu
Zae = Prae-(Piae/tand(phizae));
Pae = Prae-(Piae/tand(phipae));
numZae = [1 -Zae];
denPae = [1 -Pae];
TfZPae = tf(numZae,denPae);
Kaae = (1.1+(5.35*0.01))/norm(((polyval(numZae,Poleae)/polyval(denPae,Poleae))*(polyval(numae,Poleae)/polyval(denae,Poleae))));
figure()
plot(Prae,Piae,'p')
hold on
plot(Prae,-Piae,'p')
rlocus(Kaae*TfZPae*TFae)
figure()
margin(Kaae*TfZPae*TFae);
[GM,PM,Wp,Wg] = margin(Kaae*TfZPae*TFae);
PMae = PMtae*(180/pi)*Wg
TFfaef = feedback(Kaae*TfZPae*TFae,1);
figure()
lsim(TFfaef,ones(size(temps)),temps);
%step(TFfaef)
hold on 
plot(temps,1.02*ones(size(temps)),'-.b')
plot(temps,0.98*ones(size(temps)),'-.b')
plot(temps,1.35*ones(size(temps)),'-.b')
xline(1.25);
[numapae,denapae] = tfdata(Kaae*TfZPae*TFae,'v');
denapae = denapae(1:6);
Kaccae = abs((polyval(numapae,0)/polyval(denapae,0)));
errper = 1/Kaccae
% figure()
% test = lsim(TFfaef,temps,temps);
% plot(temps,temps-test)
% hold on
% plot(temps,Errptae*1.02*ones(size(temps)),'-.b')
% plot(temps,Errptae*0.98*ones(size(temps)),'-.b')
% xline(1.25);

%% PI
PZpiae = Prae/15;
PPpiae = 0;
numpiZae = [1 -PZpiae];
denpiPae = [1 -PPpiae];
TfZPiae = tf(numpiZae,denpiPae);
Krae = 1/abs(((polyval(numpiZae,Poleae)/polyval(denpiPae,Poleae))*(polyval(numapae,Poleae)/polyval(denapae,Poleae))));
[b,a] = cheby1(6,1,[50 60],'stop','s');
cp = tf(b,a);
TFfaepi = feedback(Kaae*TfZPae*TFae*Krae*TfZPiae,1);
figure()
margin(Kaae*TfZPae*TFae*Krae*TfZPiae) 
figure()
plot(Prae,Piae,'p')
hold on
plot(Prae,-Piae,'p')
rlocus(Kaae*TfZPae*TFae*Krae*TfZPiae)
figure()
%lsim(TFfaepi,ones(size(temps)),temps);
step(TFfaepi)
hold on 
plot(temps,1.02*ones(size(temps)),'-.b')
plot(temps,0.98*ones(size(temps)),'-.b')
plot(temps,1.3*ones(size(temps)),'-.b')
xline(1.25);
[GM,PM,Wp,Wg] = margin(Kaae*TfZPae*TFae*Krae*TfZPiae);
PMae = PMtae*(180/pi)*Wg;
%figure()
% test = lsim(TFfaepi,temps,temps);
%plot(temps,temps-test)
% endfigure()
[b,a] = cheby1(1,1,[49 61],'stop','s');
cp = tf(b,a);