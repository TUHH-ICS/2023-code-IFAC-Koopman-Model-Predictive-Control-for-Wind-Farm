function turbInputSet = controlSet_sowfa_2turb_alm_turbl(Wp)
    turbInputSet = struct();
    
    if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
end
    
    addpath([fileparts(which(mfilename)) '/LES_database']); % Add LES database
    loadedDB = load('DB_sowfa_2turb_alm_turbl.mat');
    
    turbInputSet.t = [loadedDB.turbInput.t];
    turbInputSet.phi = [loadedDB.turbInput.phi];
    turbInputSet.CT_prime = [loadedDB.turbInput.CT_prime];
    turbInputSet.interpMethod = 'lin'; % Linear interpolation over time
    
    if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
        error('Number of turbines in layout does not match your controlSet.');
    end
end