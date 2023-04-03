clear; close all; clc;
rng(0);

% Parameters from base workspace, read after setting a break point in
% Wind.m, l.678
U = 10;
Iref = 0.1400;
Lambda1 = 63; %63.0000;
T =  660000; % This was changed
t = 0 : 1:T;
%Vref = 50;

%% From Wind.m, ll 678-709

% NTM values
sigma1 = Iref*(0.75*U+5.6);
sigma2 = 0.8*sigma1;
sigma3 = 0.5*sigma1;
Lu = 8.1*Lambda1;
Lv = 2.7*Lambda1;
Lw = 0.66*Lambda1;

% Wave number vector
L = U * T;
N = length(t);
m = ifftshift(-N/2:N/2-1);
k = 2*pi*m/L;

% Spectrum
Fu = sigma1^2 * 4*Lu/U./(1+6/(2*pi)*abs(k)*Lu).^(5/3);
Fv = sigma2^2 * 4*Lv/U./(1+6/(2*pi)*abs(k)*Lv).^(5/3);
Fw = sigma3^2 * 4*Lw/U./(1+6/(2*pi)*abs(k)*Lw).^(5/3);
n = randn([3,length(k)]) + sqrt(-1)*randn([3,length(k)]);
dZ = sqrt(2*pi*[Fu;Fv;Fw]/L) .* n;

% IFFT
u = N*real(ifft(dZ(1,:)));
v = N*real(ifft(dZ(2,:)));
w = N*real(ifft(dZ(3,:)));
u = (u -mean(u)) * sigma1/std(u) + U;
v = (v -mean(v)) * sigma2/std(v);
w = (w -mean(w)) * sigma3/std(w);

%% Plot timeseries
figure(1);
plot(t,u,t,v,t,w);
legend({'u','v','w'}, 'Location', 'NorthEast')
axis tight; grid on;
ylabel('time [s]');
xlabel('Wind [m/s]');

%% Plot histogram
% alternativly, data from previous runs can be loaded
%load('testNTWind_18_1b_660s');
%load('testNTWind_18_1b_400s');
% u = Wind1VelX;
% v =  Wind1VelY;
% w = Wind1VelZ;
% t = Time;

windAbs = sqrt(u.^2 + v.^2 + w.^2)';
windAbs15 = windAbs(t>60);
figure(2);
h = histogram(windAbs15,'Normalization','pdf');
hold on
y = min(windAbs15):0.1:max(windAbs15);
mu = mean(windAbs15);
sigma = std(windAbs15);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'r','LineWidth',1.5);

pdWbl = fitdist(windAbs15,'wbl');
ypdW = pdf('wbl', windAbs15, pdWbl.A, pdWbl.B);
plot(windAbs15,ypdW,'k.','LineWidth',1.5);
axis tight; grid on;
legend({'hist(|V|)','Normal dist.','Weibull dist.'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
ylabel('Normalized occurence [-]');
xlabel('Wind [m/s]');