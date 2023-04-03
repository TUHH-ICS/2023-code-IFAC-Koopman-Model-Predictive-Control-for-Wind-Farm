clc; clear; close all; restoredefaultpath;
addpath('WFSimCode');
addpath('KoopmanIODMD');

%% 1. Create or load open loop test data with WFSim

% Inputs for WFSim_demo
R = 1e-6; % Weights on control input changed (J = sum(e'Qe + dU'Rdu)
refstairs = 0; %refstairs: Set reference to stairs
measured = 1; %measured (for Koopman model): Use measured values as feedback
KoopmanStates = 6; %KoopmanStates (for Koopman model):  Number of Koopman states
PolyLiftingFunction = 0; %PolyLiftingFunction(for Koopman model): Higher order lifting function 2,4,6,10,12,14,18,24
controller = 0; %controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
Vinf = 8;
vinfstr = '';

% Run WFSim_demo in open loop for identification set
ControlSetStr = 'sowfa_2turb_yaw_noise_step';
WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);

% RunWFsim demo for validation set
ControlSetStr = 'sowfa_2turb_yaw_steps_Ct_comb';
WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);

%Run WFSimdemo for identification for effective wind speeds
ControlSetStr = 'sowfa_2turb_alm_turbl_AllComb';
WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);

% Run WFSimdemo for open loop yaw step
ControlSetStr = 'sowfa_2turb_yaw_steps';
WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);


runWFSimAnimationPlots;

%% 2. Generate Koopman Sys Id for windfarm controller

% Inputs
yawmode = 1; %0 for ct (pitch/torque) control, 1 for additional yaw control
filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps'; %'Vinf8dot0_sowfa_2turb_alm_turbl.mat'; %
noStates = 6; %number of Koopman states 
useVal = 1; % use extra data for validation
percentTrain = .6;
NTMstr = ''; % normal turbulent wind (if empty constant free flow wind)NTMstr
PolyVec = 0;
t0 = 240;

% Run Koopman main function for WFSim simulation environent
%MainWfSim(yawmode,filenameId,filenameVal,noStates,useVal,percentTrain);
MainWfSimVinYawPowerOut(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0);
MainWfSimVinYaw(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on
yawmode = 0; %0 for ct (pitch/torque) control, 1 for additional yaw control

filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb';
filenameVal = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb';
MainWfSimVin(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on


%% 3. Evaluate quality of Koopman Sys ID in WFsim in closed loop with MPC
% Inputs

R = 1e-6; % Weights on control input changed (J = sum(e'Qe + dU'Rdu)
refstairs = -1; %refstairs: Set reference to stairs. If -1
measured = 0; %measured (for Koopman model): Use measured values as feedback
controller = 2; %controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl'; 
PolyLiftingFunction = 0; %PolyLiftingFunction(for Koopman model): Higher order lifting function 2,4,6,10,12,14,18,24
KoopmanStates = 6; %KoopmanStates (for Koopman model):  Number of Koopman states
vinfStr = '';%'_Vinf'; %_Vinf'; %alternativ '';
kTimestep =  2;
sol_update = 0;
OptYaw = 0;
KoopAIC = 1;

% Run WFSim_demo with AIC with wind speed estimates with yaw == 0
 [~,JR0,JQ0]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
    controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC);
fig = gcf;
fig.Name = sprintf('TimeSeries_OptYaw%d_KoopAIC%d',OptYaw,KoopAIC);

% Run WFSim_demo with AIC +WRC with power estimates (thrust and yawmates with yaw >= 0 for power reference
% greater greedy power
OptYaw = 1;
 [~,JR1,JQ1]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
    controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC);
fig = gcf;
fig.Name = sprintf('TimeSeries_OptYaw%d_KoopAIC%d',OptYaw,KoopAIC);

% control)
KoopAIC = 0;
 [sol_array,JR2,JQ2,fileName]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
    controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC);
fig = gcf;
fig.Name = sprintf('TimeSeries_OptYaw%d_KoopAIC%d',OptYaw,KoopAIC);

% 
% Rvec = [1,1,0.1];
% [sol_array,JR2,JQ2,fileName]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
%     controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec);
% fig = gcf;
% fig.Name = sprintf('TimeSeries_OptYaw%d_KoopAIC%d',OptYaw,KoopAIC);


%% Copy png and eps in one folder

% Manually set file names for paper plots
origDir = fullfile(pwd,'DataT2OLWFSim','Vinf8dot0_CL_K06_P0');
paperDir = fullfile(pwd,'figPaper');

% Manually set file names for paper plots
filenames = {'R1dot0e_neg06_Update0','R1dot0e_neg06_KoopAIC1_Update0',...
    'R1dot0e_neg06_KoopAIC0_Update0'};
filetype  = {'.png','.eps'}; 

for idx = 1: length(filenames)
    for idxT = 1: length(filetype)
        SOURCE = fullfile(origDir,[filenames{idx},filetype{idxT}]);
        copyfile(SOURCE, paperDir)
    end
end


% test different Rvec
% Rvec = [1,1,0.001];
% [sol_array,JR2,JQ2,fileName]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
%     controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec);
% fig = gcf;
% fig.Name = sprintf('TimeSeries_OptYaw%d_KoopAIC%d',OptYaw,KoopAIC);
% 
% Rall = [10^-8, 10^-6, 10^-4, 10^-3];
% 
% Rvec = [1,1,0.1];
% 
% for idx = 1: length(Rall)
%     R = Rall(idx);
% [sol_array,JR2(idx),JQ2(idx),fileName]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
%     controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec);
% end
