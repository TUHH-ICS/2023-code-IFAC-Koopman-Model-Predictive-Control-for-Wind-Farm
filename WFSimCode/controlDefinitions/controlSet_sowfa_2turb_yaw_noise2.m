function turbInputSet = controlSet_sowfa_2turb_yaw_noise2(Wp)
% controlSet_sowfa_2turb_yaw_steps combines controlSet_sowfa_2turb_yaw_alm_turbl
% and controlSet_sowfa_2turb_yaw_alm_uniform with added filtered white noise
% Testsignal for MIMO Koopman controller

if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
end
rng(0);


tempMat = repmat(0:5:30,4000,1);
temp = tempMat(:);

nNoise = 20;
lNoise = ceil(length(temp)/nNoise);

shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1);
noiseVec0 = sHMat(:);

tempPhi(1,:) = temp(:) + noiseVec0(1:length(temp))*0.5; 

phiRLim = 1;
for idx = 2:length(tempPhi(1,:)) % limit the yaw rate
    aDiff = tempPhi(1,idx) - tempPhi(1,idx-1);
    tempPhi(1,idx) = tempPhi(1,idx-1) + min(max(aDiff,-phiRLim),phiRLim);
end

turbInputSet.phi(1,:) = tempPhi;
turbInputSet.phi(2,:)= zeros(size(turbInputSet.phi(1,:)));


shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1)*0.3 +1.7;
noiseVec1 = sHMat(:);

shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1)*0.3 +1.7;
noiseVec2 = sHMat(:);

turbInputSet.CT_prime = max(0.2,min([noiseVec1'; noiseVec2'],2));


turbInputSet.t = 0: length(turbInputSet.CT_prime)-1;
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time



% For visualization:
% turbInputSet = controlSet_sowfa_2turb_yaw_steps;
% figure; subplot(2,1,1); plot(turbInputSet.CT_prime'); axis tight; grid on;
% subplot(2,1,2); plot(turbInputSet.phi'); axis tight; grid on;

if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
    error('Number of turbines in layout does not match your controlSet.');
end
end