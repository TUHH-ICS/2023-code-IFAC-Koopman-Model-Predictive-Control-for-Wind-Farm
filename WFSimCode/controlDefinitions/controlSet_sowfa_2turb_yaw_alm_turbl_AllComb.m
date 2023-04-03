function turbInputSet = controlSet_sowfa_2turb_yaw_alm_turbl_AllComb(Wp)
if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
    Wp.sim.h = 1;
end
    turbInputSet = controlSet_sowfa_2turb_yaw_alm_combined2;
    turbInputSet = rmfield(turbInputSet,'CT_prime');
    turbInputSet = rmfield(turbInputSet,'t');
    tmp =  turbInputSet.phi(:,60:end);
    turbInputSet = rmfield(turbInputSet,'phi');
    %loadedDB = load('CTopenloopAllComb0.01.mat');
    loadedDB = load('CToL2turbfromOLdata.mat');
    
    % nTurbs = 2;%size(loadedDB.a,1);
    len = length(loadedDB.Ct2new) * Wp.sim.h;
    turbInputSet.t = 0:Wp.sim.h:len - Wp.sim.h;
    turbInputSet.CT_prime = [loadedDB.Ct1new;loadedDB.Ct2new];
    n = ceil(len/length(tmp));
    temp2 = repmat(tmp,1,n);
    turbInputSet.phi = temp2(:,1:len);
    
end
