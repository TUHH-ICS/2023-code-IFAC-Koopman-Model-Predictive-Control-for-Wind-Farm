function turbInputSet = controlSet_sowfa_2turb_yaw_alm_combined1(Wp)

if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
end

turbInputSet1 = controlSet_sowfa_2turb_yaw_alm_turbl;
turbInputSet2 = controlSet_sowfa_2turb_yaw_alm_uniform;

k1 = 875; %te = turbInputSet1.t(k1);
k2 = 60; %t0 = turbInputSet1.t(k2);

% turbInputSet.t = [turbInputSet1.t(1:k1), turbInputSet2.t(k2:end-1) +(te +1 -t0)];
turbInputSet.phi = 3/2*[turbInputSet1.phi(:,1:k1), turbInputSet2.phi(:,k2:end)];
turbInputSet.t = 0: length(turbInputSet.phi)-1;
temp = [turbInputSet1.CT_prime(:,1:k1), turbInputSet2.CT_prime(:,k2:end)];
meanT = mean(temp,2);
% stdT = std(temp);
turbInputSet.CT_prime  = min(max(((temp-meanT)*0.1 + diag([0.2,1.4])*meanT),0.2),2) ;
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time

if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
    error('Number of turbines in layout does not match your controlSet.');
end
% figure; subplot(2,1,1); plot(turbInputSet.CT_prime'); axis tight; grid on;
% subplot(2,1,2); plot(turbInputSet.phi'); axis tight; grid on;
end