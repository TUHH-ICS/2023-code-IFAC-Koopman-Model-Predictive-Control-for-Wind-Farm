function [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val] ...
    = eDMD_RT_Vin_P(states,stateName,statesvalid,poly,n,yawmode)
%[sys_red,FITje,Xd,Xd_p,xsim,FITje_val,fig1,xo,Koop]
%[sys_red,FITje,U,S,V,X,X_p,x]
% eDMD_RTUpdate builds a state space model from the states
% and input gathered in the simulation 

if nargin < 6
    yawmode = 0;
end
 
% Size states, inputs, outputs
nx = 2; % number of original states: Ur1,Ur2 
ny = 1; % number outputs:  sum(P1,P2 (no lifted states)
nxl = n; % number of original and lifted states: Ur1,Ur2 plus lifted states
nxle = nxl + ny; % index last state in vector
nPsi = size(states,1); % number of lifted states + outputs + inputs
nu = nPsi - nxle; % number inputs: ct1, ct2,v1 

t0 = 1; % Start point; This is just for testing: 
tend = size(states,2);
tendVal = size(statesvalid,2);

% Identification set
psi_xk_1 = states(ny+1:end,t0:(tend-1)); % delayed states and inputs [x_k-1;u_k-1]
out_xk_1 = states(1:ny,t0:(tend-1)); % delayed outputs [x_k-1;u_k-1]
X_xk = states(ny+1:nxle,(t0+1):tend); % States x_k 
psi_xk = [X_xk; out_xk_1]; % [x_k+1;y_k] % States and delayed outputs 
Inputs = psi_xk_1(nxl+1:end,:); % delayed inputs [x_k-1;u_k-1]
Deterministic = out_xk_1; % outputs (ToDo Should this line be removed?)

% Validation set
psi_xk_1_val = statesvalid(ny+1:end,t0:(tendVal-1)); % delayed states and inputs [x_k-1;u_k-1]
out_xk_1_val = statesvalid(1:ny,t0:(tendVal-1)); % delayed outputs [x_k-1;u_k-1]
% X_xk_val = statesvalid(ny+1:nxle,(t0+1):tendVal); % States x_k 
% psi_xk_val = [X_xk_val; out_xk_1_val]; % [x_k+1;y_k] % States and delayed outputs 
Inputs_val = psi_xk_1_val(nxl+1:end,:); % delayed inputs [x_k-1;u_k-1]
Deterministic_val = out_xk_1_val; 

% Psi_k = K * Psi_k_1 -> K = (Psi_k*Psi_k_1^T) * pinv(Psi_k_1*Psi_k_1^T)  =
%     K = A * pinv(G)
% Psi_k'  = Psi_k_1'* K' -> K' = pinv(psi_xk_1 * psi_xk_1') * psi_xk_1 * psi_xk'
%     K' = pinv(G) * A
Koop.G = psi_xk_1 * psi_xk_1'; % Update G (GT
Koop.A = psi_xk_1 * psi_xk';   % Update A (At)
Koop.K = pinv(Koop.G) * Koop.A;% Update Koopman operator
K = Koop.K'; K_P = K((nxl+1):end,:);
K_V = K(1:nx,:);

approxA = K(1:nxl,1:nxl); % system matrix
approxB = K(1:nxl,nxl+1:end);
approxC = K(nxl+1:end,1:nxl);  %[eye(ny),zeros(ny,nx-ny)];
approxD = K(nxl+1:end,nxl+1:end); %zeros(ny,nu);
sys_red = ss(approxA,approxB,approxC,approxD,1);

% sysC = d2c(sys_red,'tustin');
% max(real(eig(sysC)))

xo = states(1:nxle,t0);% for States Ur1, Ur2 and lifted states
[ysim1,~,xsim1] = lsim(sys_red, Inputs',[],xo(ny+1:end)); %FITjeLin = vaf(Deterministic(1:end-1),ysim1(1:end-1));

structPrev.Ur1_prev1 = psi_xk_1(1,1);
structPrev.Ur2_prev1 = psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim = nan(size(psi_xk_1));
ysim = nan(ny,length(psi_xk_1)-1);
ysim(:,1) = K_P * psi_xk_1(:,1);
xsim(1:nx,1) = K_V * psi_xk_1(:,1);
xsim(:,1) = K_psi([xsim(1:nx,1);Inputs(:,1)],poly,n,structPrev,yawmode);

for idx = 1: length(psi_xk_1)-1
    ysim(:,idx) = K_P * xsim(:,idx);
    xsim(1:nx,idx+1) = K_V * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = K_psi([xsim(1:nx,idx+1); Inputs(:,idx+1)], poly,n,structPrev,yawmode);
end
FITje = vaf(Deterministic(1:end-1),ysim);

% xo_val = statesvalid(1:nx,t0);%for States Ur1,Ur2 and lifted states
% ysim_val = lsim(sys_red, Inputs_val',[],xo_val);

structPrev.Ur1_prev1 = psi_xk_1_val(1,1);
structPrev.Ur2_prev1 = psi_xk_1_val(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

xsim_val = nan(size(psi_xk_1_val));
ysim_val = nan(ny,length(psi_xk_1_val)-1);
ysim_val(:,1) =  K_P * psi_xk_1_val(:,1);
xsim_val(1:nx,1) = K_V * psi_xk_1_val(:,1);
xsim_val(:,1) = K_psi([xsim_val(1:nx,1);Inputs_val(:,1)],poly,n,structPrev,yawmode);

for idx = 1: length(xsim_val)-1
    ysim_val(:,idx) = K_P * xsim_val(:,idx);
    xsim_val(1:nx,idx+1) = K_V * xsim_val(:,idx);
    [xsim_val(:,idx+1),~,structPrev] = K_psi([xsim_val(1:nx,idx+1); Inputs_val(:,idx+1)], poly,n,structPrev,yawmode);
end

FITje_val = vaf(Deterministic_val(:,1:end-1),ysim_val);

tmp = strsplit(stateName,';');
sys_red.StateName = tmp(ny+1:nxle); % state names
sys_red.OutputName = tmp(1:ny);
%sys_red.InputName = tmp(nxle+1: nxle + nu);

if nu == 4
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1';'Vin'};
elseif nu == 3 && yawmode
    sys_red.InputName  =  {'Ct1';'Ct2';'phi1'};
elseif nu == 3
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


