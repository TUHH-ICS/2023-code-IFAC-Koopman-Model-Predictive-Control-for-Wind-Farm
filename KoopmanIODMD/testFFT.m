load('InputsInputsVal.mat','Inputs','Inputs_val')

load('Vinf8dot0_sowfa_2turb_yaw_alm_turbl_AllComb.mat')
Inputs_val= phi1;

Fs = 1;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(Inputs_val);             % Length of signal

if mod(L,2)
    L = L-1;    
end
t = (0:L-1)*T;        % Time vector


X = Inputs_val(1,1:L);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure;
plot(f(2:end),P1(2:end)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

nNoise = 10;
lNoise = ceil(L/nNoise);
rng(0);
shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1);
noiseVec1 = sHMat(:);

shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1);
noiseVec2 = sHMat(:);

turbInputSetID.CT_prime = [noiseVec1'; noiseVec2'];


X = noiseVec1(1:L);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
hold on
plot(f(2:end),P1(2:end),'r--') 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

