function [ysim,FITje] = eDMDsimulate(states,Input_u, K,poly,n,yawmode)
% eDMD_RTUpdate generates simulated output from the states 

if nargin < 5
    yawmode = 0;
end

% Size states, inputs, outputs
nx = n; % number of lifted states: Ur1,Ur2 plus lifted states
ny = 2; % number outputs:  Ur1,Ur2 (no lifted states)
psi_xk_1 = states(:,1:end-1); % States delayed
psi_xk = states(:,2:end); % States  states(1:end-1,t0:tend-1); 
Deterministic = psi_xk(1:ny,:); 
Ky = K(1:ny,:); % Get the relevant rows of the Koopman matrix

structPrev.Ur1_prev1 = psi_xk_1(1,1);
structPrev.Ur2_prev1 = psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim = nan(size(psi_xk_1));
xsim(1:ny,1) = Ky * psi_xk_1(:,1); % states: 2 + 6 + 3: P1, P2, Ur1, Ur2, n
xsim(:,1) = K_psi([xsim(1:ny,1);Input_u(:,1)],poly,n,structPrev,yawmode);
for idx = 1: length(psi_xk_1)-1
    xsim(1:ny,idx+1) = Ky * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = K_psi([xsim(1:ny,idx+1);Input_u(:,idx+1)], poly,n,structPrev,yawmode);
end
ysim = xsim(1:ny,:);
FITje = vaf(Deterministic,ysim);

