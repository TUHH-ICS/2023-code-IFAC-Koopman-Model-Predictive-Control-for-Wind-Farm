function turbInputSet = controlSet_sowfa_2turb_alm_turbl_AllComb(Wp)
    turbInputSet = struct();
    %loadedDB = load('CTopenloopAllComb0.01.mat');
    loadedDB = load('CToL2turbfromOLdata.mat');
    nTurbs = 2;%size(loadedDB.a,1);
    turbInputSet.t = 0:Wp.sim.h:length(loadedDB.Ct2new)-1;
    turbInputSet.phi = zeros(nTurbs,length(turbInputSet.t));
    turbInputSet.CT_prime = [loadedDB.Ct1new;loadedDB.Ct2new];
    turbInputSet.interpMethod = 'lin'; % Linear interpolation over time
end
