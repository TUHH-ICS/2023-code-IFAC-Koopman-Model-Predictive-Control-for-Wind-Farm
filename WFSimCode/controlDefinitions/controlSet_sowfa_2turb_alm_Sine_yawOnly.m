function turbInputSet = controlSet_sowfa_2turb_alm_Sine_yawOnly(Wp)
    turbInputSet = struct();
    
    addpath([fileparts(which(mfilename)) '/LES_database']); % Add LES database
    %loadedDB = load('DB_sowfa_2turb_alm_uniform.mat');
    
    period = 2500;
NumPeriod = 1;
wmin = 1/(40*1.5);
wmax = 1/(40);
Usin = idinput([period 2 NumPeriod],'sine',[wmin wmax],[-40,0],[5,10,1]);

    
    turbInputSet.t = 1:1:length(Usin);%[loadedDB.turbInput.t];
    turbInputSet.phi = [Usin(:,1)';Usin(:,2)';];%[loadedDB.turbInput.phi];
    turbInputSet.CT_prime = 2*ones(size(turbInputSet.phi)); %[loadedDB.turbInput.CT_prime];
    turbInputSet.interpMethod = 'lin'; % Linear interpolation over time
    
    if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
        error('Number of turbines in layout does not match your controlSet.');
    end
end

% rate_gamma1 = diff(Usin(:,1)); rate_gamma2 = diff(Usin(:,2));
%figure;
% plot(Usin(:,1));
% hold on
% plot(Usin(:,2));
% hold off

