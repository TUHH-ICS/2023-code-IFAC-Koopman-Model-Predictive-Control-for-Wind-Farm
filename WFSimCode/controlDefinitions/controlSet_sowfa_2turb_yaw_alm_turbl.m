function turbInputSet = controlSet_sowfa_2turb_yaw_alm_turbl(Wp)

if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
end

turbInputSet = struct();

addpath([fileparts(which(mfilename)) '/LES_database']); % Add LES database
loadedDB = load('DB_sowfa_2turb_yaw_alm_turbl.mat');

turbInputSet.t = [loadedDB.turbInput.t];
turbInputSet.phi = [loadedDB.turbInput.phi]*(1);
turbInputSet.phi = [0*turbInputSet.phi(:,1:31),turbInputSet.phi(:,32:end)];
turbInputSet.CT_prime = [loadedDB.turbInput.CT_prime];
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time

if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
    error('Number of turbines in layout does not match your controlSet.');
end
% figure; subplot(2,1,1); plot(turbInputSet.CT_prime'); axis tight; grid on;
%  subplot(2,1,2); plot(turbInputSet.phi'); axis tight; grid on;

end