close all
clear all
clc
% Téléscope B

%% Éléments demandés antenne B
BWbz = 10; %rad/s
PMbz = 50; % deg +-5
% Réponse à une rampe
Errrtbz = 0.005;
% CRITÈRES DE SÉCURITÉ
GMbz = 15; %dB
ATTvbz = -12; %dB
% B7 à calculer
%Critère finaux
Tsbz = 14; %s

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