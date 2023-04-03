% Code based on
% Model Predictive Control For Wake Steering: a Koopman Dynamic Mode
% Decomposition Approach
% Master Thesis Dissertation
% Author: Nassir Rodrigues Cassamo
% Supervisors: Professor Jan-Willem Van Wingerden and Professor Jo?o Sousa
% First official commit (1.0): December 2020

%% RELEVAT INFORMATION
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamic Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data

%% (0) INITIALISE
% define type of simulation
clc; clear; close all;

% directory for identification and validation data, yaw and pitch
yawmode = 0; %0 for ct (pitch/torque) control, 1 for additional yaw control

N_Prediction_Horizon = 30;  % Number of step of prediction horizon
Library_Update_Prediction_Horizon = 5;

% qLMPC options
options.nx = 2;
options.ni = 2;
options.Library_Update_Prediction_Horizon = Library_Update_Prediction_Horizon;
options.N_Prediction_Horizon = N_Prediction_Horizon;
options.MPC_enable = 0;
options.qLMPC_iterations = 1;
options.solver = 2;

filenameId = 'Vinf6dot3_sowfa_2turb_alm_turbl_AllComb.mat';%'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat';
filenameVal = 'Vinf8dot5_sowfa_2turb_alm_turbl.mat';
% filenameId = 'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat';
% filenameVal = 'Vinf8dot5_sowfa_2turb_alm_turbl.mat';
NTMstr = 'NTM'; %''

% filenameId = 'Vinf8dot5_sowfa_2turb_yaw_alm_uniform'; %'Vinf8dot5_sowfa_2turb_yaw_alm_combined.mat'; %'Vinf8dot5_sowfa_2turb_yaw_alm_turbl_AllComb.mat';
% filenameVal = 'Vinf8dot5_sowfa_2turb_yaw_alm_uniform'; %'Vinf8dot5_sowfa_2turb_yaw_alm_combined.mat';

% User input
noStates = [6,12,12,14,18,24]; %number of Koopman states [12,12,14,18,24]
useVal = 0; % use extra data for validation
percentTrain = .6;
t0 = 1;

PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
idxPolyOn = find(noStates == 12, 1, 'first');
PolyVec(idxPolyOn) = 1;

detrendingstates = 0; %1 to take mean flow and consider turbulent fluctuations
method = -2; %0: DMD ; 1:DMDc; 2:IODMD; 3:EXTIODMD -1: EDMD
koopman = 1; %to add deterministic states to flow field data

% Simulation characterisitc (resampling)
dt = 1; %time sampling (s)

KoopmanDir = fileparts(mfilename('fullpath')); % Get Koopman directory
pathFunctions = genpath(fullfile(KoopmanDir,'Functions'));

addpath(fullfile(KoopmanDir,'2.DYNAMIC_MODE_DECOMPOSITION'),...
    pathFunctions);
codedir = mfilename('fullpath');
parentdir = fileparts(fileparts(codedir));

% DEFINE MAIN output data DIRECTORY TO STORE ALL RESULTS
DataOut = fullfile(parentdir,'DataInOutWfSim');
if ~isfolder(DataOut)
    mkdir(DataOut)
end

if isempty(NTMstr)
    DataIn = fullfile(parentdir,'DataT2OLWFSim'); %,...
else
    DataIn = fullfile(parentdir,'DataT2OLWFSim','Vinf8dot0NTM_OL_Ct');
end

if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

%% Loop for  'Vinf8dot5_sowfa_2turb_yaw_alm_combined2.mat'
if yawmode == 0 %pitch control
    
    dirName = fullfile(DataIn,filenameId);
    tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    
    if useVal == 1
        dirNameVal = fullfile(DataIn, filenameVal);
        tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
            'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    end
    
else
    
    dirName = fullfile(DataIn,filenameId);
    tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    if useVal == 1
        dirNameVal = fullfile(DataIn, filenameVal);
        tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
            'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    end
end

%% (1) ASSESS DATA
% easy solution to augment u flow field data matricx with other flow
% field data
fieldCell = fieldnames(tmpId);
useVal = exist('tmpVal','var') && useVal;

for idx = 1:length(fieldCell)
    if useVal
        tmp = [tmpId.(fieldCell{idx}), tmpVal.(fieldCell{idx})]; %#ok<NASGU>
    else
        tmp = [tmpId.(fieldCell{idx})]; %#ok<NASGU>
    end
    eval([fieldCell{idx},'= tmp;']);
end

if useVal
    tident = length(tmpId.Ct1);
else
    tident = floor(length(Ur1)*percentTrain);
end
tval = tident+1;

Wp.turbine = struct(...
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747],... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.40 ... % Turbine force scaling
    );

Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );

nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx);
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);

lenUinf = length(tmpId.v);
lenUr = length(tmpId.Ur1);
len12 = lenUinf/lenUr;

utemp = reshape(u,size(u,1),len12,lenUr);
utempT = squeeze(utemp(nTx,nTy,:));
utempM = squeeze(reshape(u,size(u,1)*len12,lenUr));

vtemp = reshape(v,size(u,1),len12,lenUr);
vtempT = squeeze(vtemp(nTx,nTy,:));
vtempM = squeeze(reshape(v,size(u,1)*len12,lenUr));

if method == 0 % All measurements in QQuv
    QQ_u1 = utempM(:,t0:tident); %'; %#ok<*USENS>
    QQ_v1 = vtempM(:,t0:tident); %'
    
    valid.QQ_u1 = utempM(:,tval:end); %utempT(tval:end)';
    valid.QQ_v1 = vtempM(:,tval:end); %vtempT(tval:end)';
else % only measurements at first turbine
    QQ_u1 = utempT(t0:tident)'; %'; %#ok<*USENS>
    QQ_v1 = vtempT(t0:tident)'; %'
    
    valid.QQ_u1 = utempT(tval:end)'; %utempT(tval:end)';
    valid.QQ_v1 = vtempT(tval:end)'; %vtempT(tval:end)';
end

QQ_p1 = p(:,t0:tident);
valid.QQ_p1 = u(:,tval:end);

%% (2) DYNAMIC MODE DECOMPOSITION
states_u = QQ_u1;%fluid flow as states, identification data set
statesvalid_u = valid.QQ_u1; %fluid flow as states, validaiton data set for comparison

if yawmode == 0
    Inputs = [Ct1(t0:tident);Ct2(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end)];
else
    Inputs = [Ct1(t0:tident);Ct2(t0:tident); phi1(t0:tident)];
    Inputs_val = [Ct1(tval:end);Ct2(tval:end); phi1(tval:end)];
end

Outputs = [PT1(t0:tident);PT2(t0:tident);FT1(t0:tident);FT2(t0:tident)];
Outputs_val = [PT1(tval:end);PT2(tval:end);FT1(tval:end);FT2(tval:end)];

% ToDo check
if detrendingstates == 1
    [states_u, meansteadystate, scalingfactor] = preprocessstates(states_u); %remove meanflow or other pre processing techniques to experiment
end

strVAF = cell(length(noStates),1);
for idx = 1: length(noStates)
    
    n = noStates(idx);
    poly = PolyVec(idx);
    
    % subdirectories in dataOut
    if poly == 1 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput_poly';
    elseif poly == 0 && yawmode == 0
        dirdmdName = 'eDMDresults_UasOutput';
    elseif poly == 1
        dirdmdName = 'eDMDresults_UasOutput_MIMO_poly';
    else
        dirdmdName = 'eDMDresults_UasOutput_MIMO';
    end
    
    dirdmd = fullfile(DataOut,dirdmdName);
    if ~ isfolder(dirdmd)
        mkdir(dirdmd);
    end
    
    
    %include non linear observables - Koopman extensions to better recover non linear dynamics
    if koopman == 1
        %         states_v = QQ_v1; %fluid flow as states, identification data set
        %         states_p = QQ_p1;
        Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
        [statesUr,stateNameUr] = koopmanstateextensionWFSim(Deterministic,poly,n);
        
        structPrev.Ur1_prev1 = Deterministic(1,1);
        structPrev.Ur2_prev1 = Deterministic(2,1);
        structPrev.M1(1) = Deterministic(1,1); % Moving mean
        structPrev.M2(1) = Deterministic(2,1);
        structPrev.k = 1;
        
        Xaug2a = nan(size(statesUr));
        for idxK = 1: length(Deterministic)
            [Xaug2a(:,idxK),st1,structPrev] = koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbine = [QQ_u1; QQ_v1];
        
        if method == 0
            [statesUV,stateName] = koopmanstateextension(QQ_u1, QQ_v1);
            states = [Xaug2a; statesUV]; %(3:end,:)
        else
            [statesUV,temp] = koopmanstateextensionWFSim(windAtTurbine,12,n); %poly
            stateNameUV = regexprep(temp,{'Ur1','Ur2','M'},{'u','v','Muv'});
            states = [Xaug2a; statesUV];  %; %[statesUr]; %
            stateName = [stateNameUr,';',stateNameUV];
        end
        
        % temp = sprintf('; u%02d',1:size(states_u,1));
        % stateName = [stateNameUr,sprintf('; u%02d',1:size(states_u,1))];
        
        % statesvalid_v = valid.QQ_v1; %fluid flow as states, identification data set
        % statesvalid_p = valid.QQ_p1;
        
        Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
        windAtTurbineVal = [valid.QQ_u1; valid.QQ_v1];
        statesvalidUr  = koopmanstateextensionWFSim(Deterministic_val,poly,n);
        if method == 0
            statesUV_val = koopmanstateextension(valid.QQ_u1, valid.QQ_v1);
            statesvalid = [statesvalidUr; statesUV_val];
            
        else
  
            statesUV_val = koopmanstateextensionWFSim(windAtTurbineVal,12,n);
            statesvalid = [statesvalidUr;statesUV_val];
        end
        
        
    else
        states = states_u;
        statesvalid = statesvalid_u;
        temp = sprintf('u%d;',1: size(states,1));
        stateName = temp(1:end-1);
    end
    
    Vinf = QQ_u1(1,1);
    folderName = strrep(sprintf('Vinf%2.1f%s_states%02d',Vinf,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
    fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); %stateName);
    
    dirFig = fullfile(dirdmd,folderName); % this depends on vinf
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end
    
%     cellStName = regexp(stateName,';','split');
%     cellStName(end+1:end+3) ={'ct1','ct2','U2p'};
%     aTable = array2table([states(:,1:end-1)',Inputs(:,1:end-1)',states(1,2:end)'],...
%         'VariableNames',cellStName);
%     aMatrix = [states(:,1:end-1)',Inputs(:,1:end-1)',states(2,2:end)'];
    
    Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
    Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
    [sys_red,FITje,Xd,Xd_p,x,FITje_val,fig1,xo,Koop] = ...
        eDynamicmodedecomposition(states,statesvalid,Inputs,Outputs,...
        Deterministic,Deterministic_val,Inputs_val,Outputs_val,...
        method,dt,stateName,dirFig);
    
    save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
    
end

fid = fopen(['VAF_',strrep(filenameId,'.mat',''),'.txt'],'w');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');

for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF{idx});
end
fclose(fid);


return;

x0 = statesvalid(1:options.nx,1);
x1 = statesvalid(1:options.nx,2);
u0 = Inputs_val(:,1);
u1 = Inputs_val(:,2);

testu1 = Inputs_val(1,:);
testu2 = Inputs_val(2,:);

vecPredHor = 1:Library_Update_Prediction_Horizon+1;
xk = statesvalid(1:options.nx,vecPredHor+1); %kron(ones(Library_Update_Prediction_Horizon,1),[x1;u1]); % ToDo ADi Second vector [x0;u]
xk_ = statesvalid(1:options.nx,vecPredHor); % kron(ones(Library_Update_Prediction_Horizon,1),[x1;u1]);
Uf = Inputs_val(1:options.ni,vecPredHor);

Xk_ = koopmanstateextensionRT(xk_(1:options.nx,1),poly,n,structPrev);

Koop.v = 10;
Koop.P = 0;
Koop.error_norm_channels = 1:2;

t = length(statesvalid)-7;

%Dummy inputs
MPC.Uold    = []; %xk(end - options.ni+1:end);
MPC.ref     = [];
mu = [];
lcp2 = [];
Koop.P = 0;
%XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);

for idx = 1 : t  % simulation loop
    %MPC
    if(mod(idx,1000) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', idx, length(t)-1)
    end
    
    
    %     testU(1:2:2*N_Prediction_Horizon) = testu1( idx+1 : idx+N_Prediction_Horizon);
    %     testU(2:2:2*N_Prediction_Horizon) = testu2( idx+1 : idx+N_Prediction_Horizon);
    %     %[~,Koop,XXg,Uf,mu,xkoop,add,timing] = KqLMPC2(Koop,xk,xk_,MPC,testU,XX,mu,options);
    [Koop,XX,~,lcp,xkoop,add,timing]  = KqLMPCwt(Koop,xk,xk_,MPC,Uf,lcp2,options,poly,n,structPrev);
    timings(1:4,idx) = timing;
    
    Add(:,idx) = add;
    Xkoop(:,idx) = xkoop(:);
    
    vecPredHor = vecPredHor + 1;
    xk_ = xk; %XX(:,1:end-1); %
    %xk = [xk((nz+1):end);[x; u]];
    xk = statesvalid(1:options.nx,vecPredHor+1); % kron(ones(Library_Update_Prediction_Horizon,1),[x1;u1]); XX(:,2:end);%
    Uf = Inputs_val(1:options.ni,vecPredHor);
    
end

fid = fopen(['VAF_',strrep(filenameId,'.mat',''),'.txt'],'w');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');

figure; plot(1:t,statesvalid(1,8:end),1:t,Xkoop(1,:),'b:')
figure; plot(1:t,statesvalid(2,8:end),1:t,Xkoop(2,:),'r--')

for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF{idx});
end
fclose(fid);

%% Unused code
% Turbine and flow characteristics to be used
% rho = 1.20; %air density in [kg m^-3]
% D = 126.4; %Rotor Diameter used in simulations = 110 [m] %ToDo AD Check this in WFSim

