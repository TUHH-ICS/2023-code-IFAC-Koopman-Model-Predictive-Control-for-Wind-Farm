function MainWfSimVinYaw(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on

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
%% (0) INITIALISE
% define type of simulation
if ~nargin
    
    % directory for identification and validation data, yaw and pitch
    yawmode = 1; %0 for ct (pitch/torque) control, 1 for additional yaw control
    % User input
    noStates = 6; %[6,12,12,14,16,18,24];
    useVal = 1; % use extra data for validation
    percentTrain = .9;
    t0 = 340;
    
    NTMstr = '';
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step'; %Vinf8dot0_sowfa_2turb_yaw_noise2';
    filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps_Ct_comb';
       
    PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
    % idxPolyOn = find(noStates == 12, 1, 'first');
    % PolyVec(idxPolyOn) = 1;
end

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
    sepStr = strfind(filenameId,'_');
else
    DataIn = fullfile(parentdir,'DataT2OLWFSim','Vinf8dot0NTM_OL_Ct');
end
if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

%% Loop for  'Vinf8dot5_sowfa_2turb_yaw_alm_combined2.mat'
if yawmode == 0 %pitch control
    dirName = fullfile(DataIn, [filenameId(1:sepStr-1),'_OL_Ct'], filenameId);
    
    %     sepStr = strfind(filenameId,'_');
    %     dirNameId = fullfile(DataIn, [filenameVal(1:sepStr-1),'_OL_Ct'], filenameVal);
    tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    
    if useVal == 1
        sepStr = strfind(filenameVal,'_');
        
        filenameVal2 = strrep(filenameVal,'To','');
        dirNameVal = fullfile(DataIn, [filenameVal(1:sepStr-1),'_OL_Ct'], filenameVal2);
        tmpVal = load(dirNameVal,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
            'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    end
    
else
    
    dirName = fullfile(DataIn,[filenameId(1:sepStr-1),'_OL_Ct'],filenameId);
    tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
        'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
    if useVal == 1
        dirNameVal = fullfile(DataIn,[filenameVal(1:9),'_OL_Ct'],filenameVal);
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
tval = tident+t0;

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

kk = length(u)/Wp.mesh.Ny;

utemp = reshape(u,Wp.mesh.Nx,Wp.mesh.Ny,kk);
utempT = squeeze(utemp(nTx,nTy,:));
vtemp = reshape(v,Wp.mesh.Nx,Wp.mesh.Ny,kk);
vtempT = squeeze(vtemp(nTx,nTy,:));

QQ_u1 = utempT(t0:tident)'; %#ok<*USENS>
QQ_v1 = vtempT(t0:tident)';
QQ_p1 = p(:,t0:tident);

valid.QQ_u1 = utempT(tval:end)';
valid.QQ_v1 = vtempT(tval:end)';
valid.QQ_p1 = u(:,tval:end);

%% (2) DYNAMIC MODE DECOMPOSITION
states_u = QQ_u1;%fluid flow as states, identification data set
statesvalid_u = valid.QQ_u1; %fluid flow as states, validation data set for comparison

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
        
        structPrev.Ur1_prev1 = Deterministic(1,1);
        structPrev.Ur2_prev1 = Deterministic(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
        structPrev.M2(1) = structPrev.Ur2_prev1;
        structPrev.k = 1;
        
        statesUr = nan(n+2+yawmode,size(Deterministic,2));
        for idxK = 1: length(Deterministic)
            [statesUr(:,idxK),stateNameUr,structPrev] = K_psi([Deterministic(:,idxK);Inputs(:,idxK)],poly,n,structPrev,yawmode); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbine = [QQ_u1; QQ_v1];
        statesUV = koopmanstateextensionWFSim(windAtTurbine,poly,n);
        
        stateNameUV = regexprep(stateNameUr,{'Ur1','Ur2','M'},{'u','v','Muv'});
        states = [statesUr]; %;statesUV];
        stateName = [stateNameUr]; %,';',stateNameUV];
        
       
        structPrev.Ur1_prev1 = Deterministic_val(1,1);
        structPrev.Ur2_prev1 = Deterministic_val(2,1);
        structPrev.dUr1_prev1 = 0;
        structPrev.dUr2_prev1 = 0;
        structPrev.M1(1) = Deterministic_val(1,1); % Moving mean
        structPrev.M2(1) = Deterministic_val(2,1);
        structPrev.k = 1;
        
        statesvalidUr = nan(n+2+yawmode,size(Deterministic_val,2));
        for idxK = 1: length(Deterministic_val)
            [statesvalidUr(:,idxK),st1,structPrev] = K_psi([Deterministic_val(:,idxK);Inputs_val(:,idxK)],poly,n,structPrev); %koopmanstateextensionRT(Deterministic(:,idxK),poly,n,structPrev);
        end
        
        windAtTurbineVal = [valid.QQ_u1; valid.QQ_v1];
        statesUV_val = koopmanstateextensionWFSim(windAtTurbineVal,poly,n);
        statesvalid = [statesvalidUr];%;statesUV_val];
    else
        states = states_u;
        statesvalid = statesvalid_u;
        temp = sprintf('u%d;',1: size(states,1));
        stateName = temp(1:end-1);
    end
    
    Vinf = utemp(1,1,t0);
    Vinf2 = utemp(1,1,tident);  Vstr = '';
    if abs(Vinf -Vinf2) >0.01 && isempty(NTMstr)
        Vstr = sprintf('To%2.1f',Vinf2);
    end
    folderName = strrep(sprintf('Vin_Vinf%2.1f%s%s_states%02d',Vinf,Vstr,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
    fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); %stateName);
    
    dirFig = fullfile(dirdmd,folderName); % this depends on vinf
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end
    
    %     cellStName = regexp(stateName,';','split'); cellStName(end+1) ={'U2p'};
    %     aTable = array2table([states(:,1:end-1)',states(1,2:end)'],...
    %         'VariableNames',cellStName);
    %     aMatrix = [states(:,1:end-1)',Inputs(:,1:end-1)',states(2,2:end)'];
    
    % Sys id with ct1,ct2 as inputs
    [sys_red0,FITje0,Xd0,Xd_p0,xsim0,FITje_val0,xo_0,Koop0,...
        ysim,Inputs0,Deterministic0,ysim_val,Inputs_val0,Deterministic_val0]= eDMD_RT_Vin(states,stateName,statesvalid,poly,n,yawmode);
    strVal = 'Id.'; Vin = 0;
    plotEDMDinputsEffWind(ysim,Inputs0,Deterministic0, FITje0,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    plotEDMDinputsEffWind(ysim_val,Inputs_val0,Deterministic_val0,FITje_val0,dirFig,Vin,n);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    
    % Sys id with ct1,ct2,V1 as inputs
    statesVinf = [states;statesUV(1,:)];
    statesvalidVinf = [statesvalid; statesUV_val(1,:)];
    stateNameVinf = [stateName,';V1'];
    [sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,Koop,...
        ysim,InputsVinf,DeterministicVinf,ysim_val,Inputs_valOut,DeterministicVinf_val] = eDMD_RT_Vin(statesVinf,stateNameVinf,statesvalidVinf,poly,n,yawmode);
    strVal = 'Id.'; Vin = 1;
    plotEDMDinputsEffWind(ysim,InputsVinf,DeterministicVinf, FITje,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    strVal = 'Val.';
    plotEDMDinputsEffWind(ysim_val,Inputs_valOut,DeterministicVinf_val,FITje_val,dirFig,Vin,n,strVal);
    set(gcf,'name',strrep([strVal,fileName],'.stateName', ['Vin',num2str(Vin)]));
    
    save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
    strVAF0{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje0(1)), round(FITje_val0(1)),round(FITje0(2)), round(FITje_val0(2)));
    
end

%% Print output to file and save results of offline sys id
% print this to workspace and into file

fprintf('Three Inputs: ct1, ct2, V \nNo K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf('%s',strVAF{idx});
end

fprintf('\n\n');
fprintf('Two Inputs: ct1, ct2\nNo K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf('%s',strVAF0{idx});
end

fid = fopen(['VAF_Vin_',strrep(filenameId,'.mat',''),'.txt'],'w');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF{idx});
end

fprintf(fid,'\n\n');
fprintf(fid,'No K.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');
for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF0{idx});
end
fclose(fid);

Koffline = Koop;

return;

%% For non iterative update
% With v limit to 10
Koop.v = 10;
disp(Koop.v);
x0 = statesvalid(1:options.nx,1); % init value states Ur1, Ur2: Not
x = x0;
x_ = x0;
u0 = Inputs_val(:,1);
options.ni = size(u0,1);

testu1 = Inputs_val(1,:); % CT
testu2 = Inputs_val(2,:);
testu3 = Inputs_val(3,:);

vecPredHor = 1:Library_Update_Prediction_Horizon;
Uf = Inputs_val(1:options.ni,vecPredHor);
Uf_p = Inputs_val(1:options.ni,vecPredHor+1);
temp = [statesvalid(1:options.nx,vecPredHor+1); Uf_p]; temp(:);%
xk = temp(:);%kron(ones(Library_Update_Prediction_Horizon,1),[x0;u0]); % statesvalid(1:options.nx,vecPredHor+1); %
temp = [statesvalid(1:options.nx,vecPredHor);Uf]; % % statesvalid(1:options.nx,vecPredHor); %
xk_ = temp(:);%kron(ones(Library_Update_Prediction_Horizon,1),[x0;u0]); temp(:);

XX = kron(ones(N_Prediction_Horizon,1),xk);
U = [];
u = zeros(options.ni,1);
nx = length(x);
nu = length(u);
nz = nx + nu;
mu = zeros(4*N_Prediction_Horizon,1);
%Uf = zeros(2*N_Prediction_Horizon,1);

Xk_ = K_psi(xk_(:,1),poly,n,structPrev);
testU = Uf(:);

Koop.v = 100;
Koop.P = 1;
Koop.error_norm_channels = 1:2;

t = length(statesvalid)-N_Prediction_Horizon;

%Dummy inputs
MPC.Uold    = []; %xk(end - options.ni+1:end);
MPC.ref     = [];
mu = [];
lcp2 = [];

structPrev.Ur1_prev1 = statesvalid(1,1);
structPrev.Ur2_prev1 = statesvalid(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;
%XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);


for idx = N_Prediction_Horizon + 1 : t  % simulation loop
    %MPC
    if(mod(idx,1000) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', idx, t-1)
    end
    
    XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);
    tmp =  Inputs_val(:,idx : idx+N_Prediction_Horizon-1);
    testU = tmp(:);
    [Koop,XX,~,xkoop,add,timing] = KqLMPCwt(Koop,xk,xk_,MPC,testU,XX,mu,options,...
        poly,n,structPrev);
    
    KK{idx} = Koop.K;
    
    timings(1:4,idx) = timing;
    
    Add(:,idx) = add;
    Xkoop(:,idx) = xkoop;
    
    mu(isnan(mu)) = 0;
    %[time,xf] = ode45(@(t,x) Gyrosim(t,x,u),[0 Ts],x(:,end));
    
    x_ = x;
    x = statesvalid(1:options.nx,idx+1); %%xf(end,:)';
    u = Inputs_val(1: options.ni,idx+1);
    
    
    xk_ = xk;
    xk = [xk((nz+1):end);[x; u]]; %ToDo ADi 10 replaced by length([x; u]) +1
    % xsim = [xsim xf(end,:)'];
    % U = [U u];
    
    % xk = statesvalid(1:options.nx,vecPredHor+1); % kron(ones(Library_Update_Prediction_Horizon,1),[x1;u1]); XX(:,2:end);%
    Uf = Inputs_val(1:options.ni,vecPredHor);
    
end

vecT = 60:t;
figure;
cl = lines(2);
subplot(3,1,1);
plot(vecT,testu1(vecT),vecT,testu2(vecT)); ylabel('C_T[-]'); legend('WT1','WT2')
axis tight; grid on;
subplot(3,1,2)
plot(vecT,statesvalid(1,vecT),vecT,Xkoop(1,vecT),'b:'); ylabel('Ur1 [m/s]');legend('real','est')
axis tight; grid on;
subplot(3,1,3);
plot(vecT,statesvalid(2,vecT),'color',cl(2,:)); ylabel('Ur2 [m/s]');
axis tight; grid on;
hold on; plot(vecT,Xkoop(2,vecT),'r--'); legend('real','est')

figure;
xlabelCell = {'eHinf','Update','MaxRealEig','MaxRealEigUsed'};
noPl = size(Add,1);
for idx = 1:noPl
    subplot(noPl,1,idx);
    plot(Add(idx,:)); axis tight, grid on;
    ylabel(xlabelCell{idx});
end

if max(Add(1,:)) > 10
    figure;
    xlabelCell = {'eHinf','Update','MaxRealEig'};
    for idx = 1:noPl
        subplot(noPl,1,idx);
        plot(Add(idx,6000:end)); axis tight, grid on;
        ylabel(xlabelCell{idx});
    end
end

%% For iterative update
% With v limit
Koop_v = 3;
disp(Koop_v);
x0 = statesvalid(1:options.nx,1); % init value states Ur1, Ur2: Not
x = x0;
x_ = x0;
u0 = Inputs_val(:,1);
options.ni = size(u0,1);

testu1 = Inputs_val(1,:); % CT
testu2 = Inputs_val(2,:);
testu3 = Inputs_val(3,:);

vecPredHor = 1:Library_Update_Prediction_Horizon;
Uf = Inputs_val(1:options.ni,vecPredHor);
Uf_p = Inputs_val(1:options.ni,vecPredHor+1);
temp = [statesvalid(1:options.nx,vecPredHor+1); Uf_p]; temp(:);%
xk = temp(:);%kron(ones(Library_Update_Prediction_Horizon,1),[x0;u0]); % statesvalid(1:options.nx,vecPredHor+1); %
temp = [statesvalid(1:options.nx,vecPredHor);Uf]; % % statesvalid(1:options.nx,vecPredHor); %
xk_ = temp(:);%kron(ones(Library_Update_Prediction_Horizon,1),[x0;u0]); temp(:);

XX = kron(ones(N_Prediction_Horizon,1),xk);
U = [];
u = zeros(options.ni,1);
nx = length(x);
nu = length(u);
nz = nx + nu;
mu = zeros(4*N_Prediction_Horizon,1);
%Uf = zeros(2*N_Prediction_Horizon,1);

Xk_ = K_psi(xk_(:,1),poly,n,structPrev);
testU = Uf(:);

Koop = Koffline;
Koop.P = 1;
Koop.error_norm_channels = 1:2;
Koop.v = Koop_v;

t = length(statesvalid)-N_Prediction_Horizon;

%Dummy inputs
MPC.Uold    = []; %xk(end - options.ni+1:end);
MPC.ref     = [];
mu = [];
lcp2 = [];

structPrev.Ur1_prev1 = statesvalid(1,1);
structPrev.Ur2_prev1 = statesvalid(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;
%XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);
Add = nan(3,t);

for idx = N_Prediction_Horizon + 1 : t  % simulation loop
    %MPC
    if(mod(idx,1000) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', idx, t-1)
    end
    
    XX = kron(ones(N_Prediction_Horizon,1),[x;x-x_]);
    tmp =  Inputs_val(:,idx : idx+N_Prediction_Horizon-1);
    testU = tmp(:);
    [Koop,XX,~,xkoop,add,timing] = KqLMPCwt(Koop,xk,xk_,MPC,testU,XX,mu,options,...
        poly,n,structPrev);
    
    KK{idx} = Koop.K;
    
    timings(1:4,idx) = timing;
    
    Add(:,idx) = add;
    Xkoop(:,idx) = xkoop;
    
    mu(isnan(mu)) = 0;
    %[time,xf] = ode45(@(t,x) Gyrosim(t,x,u),[0 Ts],x(:,end));
    
    x_ = x;
    x = statesvalid(1:options.nx,idx+1); %%xf(end,:)';
    u = Inputs_val(1: options.ni,idx+1);
    
    xk_ = xk;
    xk = [xk((nz+1):end);[x; u]]; %ToDo ADi 10 replaced by length([x; u]) +1
    % xsim = [xsim xf(end,:)'];
    % U = [U u];
    
    % xk = statesvalid(1:options.nx,vecPredHor+1); % kron(ones(Library_Update_Prediction_Horizon,1),[x1;u1]); XX(:,2:end);%
    Uf = Inputs_val(1:options.ni,vecPredHor);
    
    % Add the real eigenvalue of
    % K = Koop.K';  % Update the truncated Koopman operator
    % ny = 2;
    %     dt = 1;
    %     approxA = K(1:noStates,1:noStates); % system matrix
    %     approxB = K(1:noStates,noStates+1:end);
    %     approxC = [eye(ny),zeros(ny,noStates-ny)];
    %     approxD = zeros(ny,nu);
    %     sysC = d2c(ss(approxA,approxB,approxC,approxD,dt),'tustin');
    %     add(3,idx) = max(real(eig(sysC)));
    
end

%% REsult of online update
vecT = 60:t;
figure;
cl = lines(2);
subplot(3,1,1);
plot(vecT,testu1(vecT),vecT,testu2(vecT)); ylabel('C_T[-]'); legend('WT1','WT2')
axis tight; grid on;
subplot(3,1,2)
plot(vecT,statesvalid(1,vecT),vecT,Xkoop(1,vecT),'b:'); ylabel('Ur1 [m/s]');legend('real','est')
axis tight; grid on;
subplot(3,1,3);
plot(vecT,statesvalid(2,vecT),'color',cl(2,:)); ylabel('Ur2 [m/s]');
axis tight; grid on;
hold on; plot(vecT,Xkoop(2,vecT),'r--'); legend('real','est')

figure;
xlabelCell = {'eHinf','Update','MaxRealEig','MaxRealEig0'};
noPl = size(Add,1);
for idx = 1:noPl
    subplot(noPl,1,idx);
    plot(Add(idx,:)); axis tight, grid on;
    ylabel(xlabelCell{idx});
end


%% Offline SysId based on data
[sys_red,FITje,Xd,Xd_p,xsim,FITje_val,xo,KoopEnd,...
    ysim,Inputs,Deterministic,ysim_val,Inputs_val,Deterministic_val] = eDMD_RT_Vin(statesvalidVinf,stateNameVinf,statesVinf,poly,n);
plotEDMDinputsEffWind(ysim,Inputs,Deterministic, FITje,dirFig,Vin,n,strVal);
plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,n);

%%
ny = 2;
K = KoopEnd.K'; Ky = K(1:ny,:);
t0 = 1;
psi_xk_1 = statesvalidVinf(:,t0:end-1); % States delayed
psi_xk = statesvalidVinf(:,t0+1:end); % States  states(1:end-1,t0:tend-1); % %X_k = states(1:nx,t0:tend-1);


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
xsim(1:ny,1) = Ky * psi_xk_1(:,1);
xsim(:,1) = K_psi([xsim(1:ny,1);psi_xk(n+1:end,1)],poly,n,structPrev);
for idx = 1: length(psi_xk_1)-1
    xsim(1:ny,idx+1) = Ky * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = K_psi([xsim(1:ny,idx+1);psi_xk_1(n+1:end,idx+1)], poly,n,structPrev);
end
ysim = xsim(1:ny,:);

FITje = vaf(statesvalidVinf(1:ny,1:end-1),ysim);

%%
ny = 2;
K = KoopEnd.K'; Ky = K(1:ny,:);
t0 = 1;
psi_xk_1 = statesvalidVinf(:,t0:end-1); % States delayed
psi_xk = statesvalidVinf(:,t0+1:end); % States  states(1:end-1,t0:tend-1); % %X_k = states(1:nx,t0:tend-1);

xo = states(1:nx,t0);% for States Ur1, Ur2 and lifted states
% [ysim,~,xsim] = lsim(sys_red, Inputs',[],xo);

structPrev.Ur1_prev1 = psi_xk_1(1,1);
structPrev.Ur2_prev1 = psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

KK(1: N_Prediction_Horizon) = KK(N_Prediction_Horizon+1);

tSim = min(length(KK),length(psi_xk_1));
xsim1 = nan(size(psi_xk_1,1),tSim);
K = KK{1}'; Ky = K(1:ny,:);
xsim1(1:ny,1) = Ky * psi_xk_1(:,1);
xsim1(:,1) = K_psi([xsim1(1:ny,1);psi_xk(n+1:end,1)],poly,n,structPrev);
for idx = 1: tSim-1
    K = KK{idx}'; Ky = K(1:ny,:);
    xsim1(1:ny,idx+1) = Ky * xsim1(:,idx);
    [xsim1(:,idx+1),~,structPrev] = K_psi([xsim1(1:ny,idx+1);psi_xk_1(n+1:end,idx+1)], poly,n,structPrev);
end
ysim = xsim1(1:ny,:);

FITje = vaf(statesvalidVinf(1:ny,1:5500),ysim(:,1:5500));
%plotEDMDinputsEffWind(ysim,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,n);


%% test with KoopK5
%load('Koop5KSim','Koop5K');
load('Koop5KSim','Koop5K','nSample','tSim');

Koop5K(tSim-nSample+2:tSim+1) = Koop5K(tSim-nSample+1);
tmpCell = Koop5K;
tmpCell(nSample:tSim) = Koop5K(1:tSim-nSample+1);
tmpCell(1:nSample-1) = tmpCell(nSample);
tmpCell(tSim: length(xsim_valF)) = tmpCell(tSim);

xsim_valF = nan(size(statesvalid));
Khat = tmpCell{1};
Ky = Khat(1:2,:);

xsim_valF(1:ny,1) = Ky * statesvalid(:,1);
xsim_valF(:,1) = K_psi([xsim_valF(1:ny,1);statesvalid(nx+1:end,1)],poly,n,structPrev);
for idx = 1: length(xsim_valF)-1
    Khat = tmpCell{idx+1};
    Ky = Khat(1:2,:);
    xsim_valF(1:ny,idx+1) = Ky * xsim_valF(:,idx);
    [xsim_valF(:,idx+1),~,structPrev]  = K_psi([xsim_valF(1:ny,idx+1);statesvalid(nx+1:end,idx+1)],poly,n,structPrev);
end
ysim_valF = xsim_valF(1:ny,2:end);

FITje_valF = vaf(Deterministic_val(:,1:length(ysim_valF)),ysim_valF);


strVal = 'OnlineId'; displayOff = 1;
plotEDMDinputsEffWind(ysim_valF,Inputs_val,Deterministic_val,FITje_valF,dirFigTest1,Vin,n,strVal,displayOff);
plotEDMDinputsEffWindUpdate(ysim_valF,Inputs_val,Deterministic_val,FITje_valF,dirFigTest1,Vin,nx,strVal,P5vec,displayOff);

%% plot sigma plots
%idxDiff = find(diff(Add(2,:)==1));
vecK = [N_Prediction_Horizon+1, length(KK)];
sysC = cell(2,1);

for idxK = 1: length(vecK)
    K = KK{vecK(idxK)}'; %';  % Update the truncated Koopman operator
    ny = 2;
    nu = 3;
    noStates = n;
    dt = 1;
    approxA = K(1:noStates,1:noStates); % system matrix
    approxB = K(1:noStates,noStates+1:end);
    approxC = [eye(ny),zeros(ny,noStates-ny)];
    approxD = zeros(ny,nu);
    sysC{idxK} = d2c(ss(approxA,approxB,approxC,approxD,dt),'tustin');
    
end

sysC{end+1} = d2c(sys_red,'tustin');
% [SV7] = sigma(sysC7); [SV79] = sigma(sysC79);

ltCell = {'-','--','-.'};
figure;
for idx = 1: length(sysC)
    sigma(sysC{idx})
    hold on;
end

hold off;

legend('sysOnline0','sysOnlineEnd', 'sysOffline');


%% Unused code
% Turbine and flow characteristics to be used
% rho = 1.20; %air density in [kg m^-3]
% D = 126.4; %Rotor Diameter used in simulations = 110 [m] %ToDo AD Check this in WFSim

