function turbInputSet = controlSet_sowfa_2turb_yaw_steps(Wp)
% controlSet_sowfa_2turb_yaw_steps combines controlSet_sowfa_2turb_yaw_alm_turbl
% and controlSet_sowfa_2turb_yaw_alm_uniform with added filtered white noise
%
% For visualization:
% turbInputSet = controlSet_sowfa_2turb_yaw_steps;
% figure; subplot(2,1,1); plot(turbInputSet.CT_prime'); axis tight; grid on;
% subplot(2,1,2); plot(turbInputSet.phi'); axis tight; grid on;

if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
    Wp.sim.h = 1;
end

temp = repmat(0:5:40,500,1);
tempPhi(1,:) = temp(:); %[zeros(1000,1); 5*ones(900,1);10*ones(900,1);...
%15*ones(900,1);20*ones(900,1); 25*ones(900,1)];

phiRLim = 1;

for idx = 2:length(tempPhi(1,:)) % limit the yaw rate
    aDiff = tempPhi(1,idx) - tempPhi(1,idx-1);
    tempPhi(1,idx) = tempPhi(1,idx-1) + min(max(aDiff,-phiRLim),phiRLim);
end

turbInputSet.phi(1,:) = tempPhi;
turbInputSet.phi(2,:)= zeros(size(turbInputSet.phi(1,:)));
turbInputSet.CT_prime = 2 *ones(size(turbInputSet.phi));
turbInputSet.t = 0: length(turbInputSet.CT_prime)-1;
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time

end
