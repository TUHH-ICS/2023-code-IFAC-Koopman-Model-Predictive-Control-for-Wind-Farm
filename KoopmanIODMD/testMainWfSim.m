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
Library_Update_Prediction_Horizon = 5; %5;
% qLMPC options
options.nx = 2;
options.ni = 2;
nz = options.nx +options.ni;
options.Library_Update_Prediction_Horizon = Library_Update_Prediction_Horizon;
options.N_Prediction_Horizon = N_Prediction_Horizon;
options.MPC_enable = 0;
options.qLMPC_iterations = 1;
options.solver = 2;


% MPC options
MPC.Q_      = blkdiag(kron(eye(N_Prediction_Horizon-1),...
    diag([5 1500 1000 .01 5 10 1 0 0 0 0 0 0 0])),...
    10*diag([5 100 100 .01 5 10 1 0 0 0 0 0 0 0]));
MPC.R_      = kron(eye(N_Prediction_Horizon),10000*eye(options.ni));
% Paper Values
MPC.Q_      = blkdiag(kron(eye(N_Prediction_Horizon-1),diag([1 120 120 .01 5 2 2 0 0 0 0 0 0 0])),10*diag([5 100 100 .01 5 10 1 0 0 0 0 0 0 0]));
MPC.R_      = kron(eye(N_Prediction_Horizon),diag([3000,750]));
MPC.C_delta = tril(kron(ones(N_Prediction_Horizon),eye(options.ni)));
MPC.Ulim    = kron(ones(N_Prediction_Horizon,1),[0.5;2]);
MPC.Caug    = kron(eye(N_Prediction_Horizon),[1 0 0 0 0 0 0 0 0 0 0 0 0 0]);
MPC.ylim    = pi/2;


% filenameId = 'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat';
% %filenameId = 'Copy_of_Vinf6dot3_sowfa_2turb_alm_turbl_AllComb.mat'; %'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat';
%
% filenameVal = 'Vinf8dot5_sowfa_2turb_alm_turbl.mat';

filenameId = 'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat'; NTMstr = '';
%filenameId = 'Copy_of_Vinf6dot3_sowfa_2turb_alm_turbl_AllComb.mat'; NTMstr = 'NTM';%'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat';
% filenameVal = 'Vinf8dot5_sowfa_2turb_alm_turbl.mat';
% filenameId = 'Vinf8dot5_sowfa_2turb_alm_turbl_AllComb.mat';
% filenameVal = 'Vinf8dot5_sowfa_2turb_alm_turbl.mat';
NTMstr = '';

% filenameId = 'Vinf8dot5_sowfa_2turb_yaw_alm_uniform'; %'Vinf8dot5_sowfa_2turb_yaw_alm_combined.mat'; %'Vinf8dot5_sowfa_2turb_yaw_alm_turbl_AllComb.mat';
% filenameVal = 'Vinf8dot5_sowfa_2turb_yaw_alm_uniform'; %'Vinf8dot5_sowfa_2turb_yaw_alm_combined.mat';

% User input
noStates = 12; %[6,12,12,14,18,24]
useVal = 1; % use extra data for validation
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

Nx = Wp.mesh.Nx;
Ny = Wp.mesh.Ny;
nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Nx);
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Ny);
lenUr = length(u)/Ny;

utemp = reshape(u,Nx,Ny,lenUr);
utempT = squeeze(utemp(nTx,nTy,:));
vtemp = reshape(v,Nx,Ny,lenUr);
vtempT = squeeze(vtemp(nTx,nTy,:));

QQ_u1 = utempT(t0:tident)'; %#ok<*USENS>
QQ_v1 = vtempT(t0:tident)';
QQ_p1 = p(:,t0:tident);

valid.QQ_u1 = utempT(tval:end)';
valid.QQ_v1 = vtempT(tval:end)';
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
    
    Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
    Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
    %include non linear observables - Koopman extensions to better recover non linear dynamics
    if koopman == 1
        %         states_v = QQ_v1; %fluid flow as states, identification data set
        %         states_p = QQ_p1;
        
        %[statesUr,stateNameUr] = koopmanstateextensionWFSim(Deterministic,poly,n);
        
        structPrev.Ur1_prev1 = Deterministic(1,1);
        structPrev.Ur2_prev1 = Deterministic(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = Deterministic(1,1); % Moving mean
        structPrev.M2(1) = Deterministic(2,1);
        structPrev.k = 1;
        
        statesUr = nan(n+2,size(Deterministic,2));
        for idxK = 1: length(Deterministic)
            [statesUr(:,idxK),stateNameUr,structPrev] = K_psi([Deterministic(:,idxK);Inputs(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbine = [QQ_u1; QQ_v1];
        statesUV = koopmanstateextensionWFSim(windAtTurbine,poly,n);
        stateNameUV = regexprep(stateNameUr,{'Ur1','Ur2','M'},{'u','v','Muv'});
        states = statesUr; %[statesUr]; %;statesUV];
        stateName = [stateNameUr]; %,';',stateNameUV];

        structPrev.Ur1_prev1 = Deterministic_val(1,1);
        structPrev.Ur2_prev1 = Deterministic_val(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = Deterministic_val(1,1); % Moving mean
        structPrev.M2(1) = Deterministic_val(2,1);
        structPrev.k = 1;
        statesvalidUr = nan(n+2,size(Deterministic_val,2));
        for idxK = 1: length(Deterministic_val)
            [statesvalidUr(:,idxK),st1,structPrev] = K_psi([Deterministic_val(:,idxK);Inputs_val(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        windAtTurbineVal = [valid.QQ_u1; valid.QQ_v1];
        statesUV_val = koopmanstateextensionWFSim(windAtTurbineVal,poly,n);
        statesvalid = [statesvalidUr]; %;statesUV_val];
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
    
    %     cellStName = regexp(stateName,';','split'); cellStName(end+1) ={'U2p'};
    %     aTable = array2table([states(:,1:end-1)',states(1,2:end)'],...
    %         'VariableNames',cellStName);
    %     aMatrix = [states(:,1:end-1)',Inputs(:,1:end-1)',states(2,2:end)'];
    
    %     [sys_red,FITje,Xd,Xd_p,x,FITje_val,fig1,xo,Koop] = ...
    %         eDynamicmodedecomposition(states,statesvalid,Inputs,Outputs,...
    %         Deterministic,Deterministic_val,Inputs_val,Outputs_val,...
    %         method,dt,stateName,dirFig);
    
    [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,fig1,xo,Koop] = eDMD_RTUpdate(states,dt,stateName,statesvalid,dirFig);
    
    
    % save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
    
end


%% For iterative update
x0 = statesvalid(1:options.nx,1); % init value states Ur1, Ur2: Not
x = x0;
x_ = x0;
u0 = Inputs_val(:,1);

testu1 = Inputs_val(1,:); % CT
testu2 = Inputs_val(2,:);
Koop.vLearn = 10;
Koffline = Koop.K;

vecPredHor = 1:Library_Update_Prediction_Horizon;
%xk = kron(ones(Library_Update_Prediction_Horizon,1),[x0;u0]); % ToDo ADi Second vector [x0;u] statesvalid(1:options.nx,vecPredHor+1); %
tmp = [statesvalid(1:options.nx,vecPredHor+1); testu1(vecPredHor+1); testu2((vecPredHor+1))];
xk = tmp(:);
% xk_ = kron(ones(Library_Update_Prediction_Horizon,1),[x0;u0]); % statesvalid(1:options.nx,vecPredHor); %
tmp = [statesvalid(1:options.nx,vecPredHor); testu1(vecPredHor); testu2((vecPredHor))];
xk_ = tmp(:);
Uf = Inputs_val(1:options.ni,vecPredHor);

XX = kron(ones(N_Prediction_Horizon,1),xk);
U = [];
u = zeros(2,1);
nx = length(x);
nu = length(u);
nz = nx + nu;
mu = zeros(4*N_Prediction_Horizon,1);
Uf = zeros(2*N_Prediction_Horizon,1);

Xk_ = K_psi(xk_(:,1),poly,n,structPrev);
testU = Uf;

Koop.v = 10;
Koop.P = 0;
Koop.error_norm_channels = 1:2;

t = length(statesvalid) - N_Prediction_Horizon;

%Dummy inputs
MPC.Uold    = []; %xk(end - options.ni+1:end);
MPC.ref     = [];
mu = [];
lcp2 = [];
Koop.P = 0;
%XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);

K0 = Koop.K';

% Koop.K = zeros(8);
% Koop.G = zeros(8);
% Koop.A = zeros(8);
% Koop.P = 0; %Koop2 = Koop;

for idx = Library_Update_Prediction_Horizon +1 : t  % simulation loop
    %MPC
    if(mod(idx,1000) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', idx, length(t)-1)
    end
    
    %
    %     testU(1:2:2*N_Prediction_Horizon) = testu1( idx+1 : idx+N_Prediction_Horizon);
    %     testU(2:2:2*N_Prediction_Horizon) = testu2( idx+1 : idx+N_Prediction_Horizon);
    %     %     %[~,Koop,XXg,Uf,mu,xkoop,add,timing] = KqLMPC2(Koop,xk,xk_,MPC,testU,XX,mu,options);
    %     [Koop,XX,~,lcp,xkoop,add,timing]  = KqLMPCwt(Koop,xk,xk_,MPC,testU,lcp2,options,poly,n,structPrev);
    %
    %Koop.v = Koop.vLearn;
    XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);
    testU(1:2:end) = testu1( idx+1 : idx+N_Prediction_Horizon);
    testU(2:2:end) = testu2( idx+1 : idx+N_Prediction_Horizon);
    [~,Koop,XXg,Uf,mu,xkoop,add,timing] = KqLMPC2(Koop,xk,xk_,MPC,testU,XX,mu,options);
    
    KK{idx} = Koop.K;
    
    timings(1:4,idx) = timing;
    
    Add(:,idx) = add;
    Xkoop(:,idx) = xkoop;
    
    mu(isnan(mu)) = 0;
    %[time,xf] = ode45(@(t,x) Gyrosim(t,x,u),[0 Ts],x(:,end));
    
    x_ = x;
    x = statesvalid(1:options.nx,idx+1); %%xf(end,:)';
    u = statesvalid(end-(options.ni-1):end,idx+1);
    
    %  x_ = x;
    %  x = xf(end,:)';
    
    xk_ = xk;
    xk = [xk((nz+1):end);[x; u]]; %ToDo ADi 10 replaced by length([x; u]) +1
    
    %     xsim = [xsim xf(end,:)'];
    %     U = [U u];
    
    % xk = statesvalid(1:options.nx,vecPredHor+1); % kron(ones(Library_Update_Prediction_Horizon,1),[x1;u1]); XX(:,2:end);%
    Uf = Inputs_val(1:options.ni,vecPredHor);
    
    % Add the real eigenvalue of
    K = Koop.K';  % Update the truncated Koopman operator
    
    %
    ny = 2;
    dt = 1;
    noStates = 6;
    approxA = K(1:noStates,1:noStates); % system matrix
    approxB = K(1:noStates,noStates+1:end);
    approxC = [eye(ny),zeros(ny,noStates-ny)];
    approxD = zeros(ny,nu);
    sysC = d2c(ss(approxA,approxB,approxC,approxD,dt),'tustin');
    add(3,idx) = max(real(eig(sysC)));
    
end


fid = fopen(['VAF_',strrep(filenameId,'.mat',''),'.txt'],'w');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');

vecT = 6:t;
figure;
cl = lines(2);
subplot(3,1,1);
plot(vecT,testu1(vecT),vecT,testu2(vecT)); ylabel('C_T[-]'); legend('WT1','WT2')
subplot(3,1,2)
plot(vecT,statesvalid(1,vecT),vecT,Xkoop(1,vecT),'b:'); ylabel('Ur1 [m/s]');legend('real','est')
subplot(3,1,3);
plot(vecT,statesvalid(2,vecT),'color',cl(2,:)); ylabel('Ur2 [m/s]');
hold on; plot(vecT,Xkoop(2,vecT),'r--'); legend('real','est')

for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF{idx});
end
fclose(fid);

figure;
xlabelCell = {'eHinf','Update','MaxRealEig'};
noPl = size(Add,1);
for idx = 1:noPl
    subplot(noPl,1,idx);
    plot(Add(idx,:)); axis tight, grid on;
    ylabel(xlabelCell{idx});
end


figure;
xlabelCell = {'eHinf','Update','MaxRealEig'};
for idx = 1:noPl
    subplot(noPl,1,idx);
    plot(Add(idx,6000:end)); axis tight, grid on;
    ylabel(xlabelCell{idx});
end

%% Unused code
% Turbine and flow characteristics to be used
% rho = 1.20; %air density in [kg m^-3]
% D = 126.4; %Rotor Diameter used in simulations = 110 [m] %ToDo AD Check this in WFSim

