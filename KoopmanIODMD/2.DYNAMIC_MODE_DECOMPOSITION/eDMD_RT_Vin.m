function [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val] ...
    = eDMD_RT_Vin(states,stateName,statesvalid,poly,n,yawmode)
%[sys_red,FITje,Xd,Xd_p,xsim,FITje_val,fig1,xo,Koop]
%[sys_red,FITje,U,S,V,X,X_p,x]
% eDMD_RTUpdate builds a state space model from the states
% and input gathered in the simulation 

if nargin < 6
    yawmode = 0;
end

% Size states, inputs, outputs
nPsi = size(states,1); % number of lifted states + inputs
nx = n; % number of lifted states: Ur1,Ur2 plus lifted states
nu = nPsi - n; %size(inp,1); % number inputs: ct1, ct2,v1 
ny = 2; % number outputs:  Ur1,Ur2 (no lifted states)

t0 = 1;
tend = round(size(states,2));
tendVal = size(statesvalid,2);

psi_xk_1 = states(:,t0:tend-1); % States delayed
psi_xk = states(:,t0+1:tend); % States  states(1:end-1,t0:tend-1); % %X_k = states(1:nx,t0:tend-1);

Inputs = psi_xk(end-(nu-1): end,:);
Deterministic = psi_xk(1:ny,:); 

Inputs_val = statesvalid((nPsi -(nu-1): nPsi), (t0:tendVal-1));
Deterministic_val = statesvalid(1:ny, (t0+1) : tendVal); 

% Psi_k = K * Psi_k_1 -> K = (Psi_k*Psi_k_1^T) * pinv(Psi_k_1*Psi_k_1^T)  =
%     K = A * pinv(G)
% Psi_k'  = Psi_k_1'* K' -> K' = pinv(psi_xk_1 * psi_xk_1') * psi_xk_1 * psi_xk'
%     K' = pinv(G) * A
Koop.G = psi_xk_1 * psi_xk_1'; % Update G (GT
Koop.A = psi_xk_1 * psi_xk';   % Update A (At)
Koop.K = pinv(Koop.G) * Koop.A;% Update Koopman operator
K = Koop.K'; Ky = K(1:ny,:);

approxA = K(1:nx,1:nx); % system matrix
approxB = K(1:nx,nx+1:end);
approxC = [eye(ny),zeros(ny,nx-ny)];
approxD = zeros(ny,nu);
sys_red = ss(approxA,approxB,approxC,approxD,1);

% sysC = d2c(sys_red,'tustin');
% max(real(eig(sysC)))

xo = states(1:nx,t0);% for States Ur1, Ur2 and lifted states
% [ysim,~,xsim] = lsim(sys_red, Inputs',[],xo);

structPrev.Ur1_prev1 = psi_xk_1(1,1);
structPrev.Ur2_prev1 = psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim = nan(size(psi_xk_1));
xsim(1:ny,1) = Ky * psi_xk_1(:,1); % states: 2 + 6 + 3: P1, P2, Ur1, Ur2, n
xsim(:,1) = K_psi([xsim(1:ny,1);psi_xk(nx+1:end,1)],poly,n,structPrev,yawmode);
for idx = 1: length(psi_xk_1)-1
    xsim(1:ny,idx+1) = Ky * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = K_psi([xsim(1:ny,idx+1);psi_xk_1(nx+1:end,idx+1)], poly,n,structPrev,yawmode);
end
ysim = xsim(1:ny,:);

FITje = vaf(Deterministic,ysim);

% xo_val = statesvalid(1:nx,t0);%for States Ur1,Ur2 and lifted states
% ysim_val = lsim(sys_red, Inputs_val',[],xo_val);
structPrev.Ur1_prev1 = statesvalid(1,1);
structPrev.Ur2_prev1 = statesvalid(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim_val = nan(size(statesvalid));
xsim_val(1:ny,1) = Ky * statesvalid(:,1);
xsim_val(:,1) = K_psi([xsim_val(1:ny,1);statesvalid(nx+1:end,1)],poly,n,structPrev);
for idx = 1: length(xsim_val)-1
    xsim_val(1:ny,idx+1) = Ky * xsim_val(:,idx);
    [xsim_val(:,idx+1),~,structPrev]  = K_psi([xsim_val(1:ny,idx+1);statesvalid(nx+1:end,idx+1)],poly,n,structPrev,yawmode);
end
ysim_val = xsim_val(1:ny,2:end);

FITje_val = vaf(Deterministic_val,ysim_val);

tmp = strsplit(stateName,';');
sys_red.StateName = tmp(1:nx); % state names
sys_red.OutputName = {'Ur1';'Ur2'};

if max(size(sys_red)) == 4
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1';'Vin'};
elseif max(size(sys_red)) == 3 && yawmode
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1'};
elseif max(size(sys_red)) == 3
    sys_red.InputName  =  {'Ct1';'Ct2';'Vin'};
else
    sys_red.InputName  =  {'Ct1';'Ct2'};
end

% Ouput names should be consistent
Xd = psi_xk_1; Xd_p = psi_xk;

% Plot for validation data
% Vin = 1;
% fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val, FITje_val,dirFig,Vin,nx);
% strVal = 'Id.';
% plotEDMDinputsEffWind(ysim,Inputs,Deterministic, FITje,dirFig,Vin,nx,strVal);


% Original code
%     Psi = [X; inp]; % Output Ur1, Ur2
%     Koop.G = Psi * Psi';
%     Koop.A = X_p * Psi';
%     Koop.K = Koop.A * pinv(Koop.G); %X_p * pinv([X;inp]);%for States Ur1,Ur2 and lifted states
% Code for debugging
% figure; 
% 
% subplot(2,1,1)
% plot(temp(1,:))
% hold on; plot(psi_xk(1,:),'--')
% subplot(2,1,1)
% plot(temp(2,:))
% hold on; plot(psi_xk(2,:),'--')
% 

% figure; 
% subplot(2,1,1)
% plot(tmp2(1,1:end))
% hold on; plot(psi_xk(1,1:end),'--'); plot(psi_xk(end,1:end),'k:'); hold off;
% subplot(2,1,2)
% plot(tmp2(2,1:end));hold on; plot(psi_xk(2,1:end),'--'); hold off;
% temp = K * psi_xk_1;


