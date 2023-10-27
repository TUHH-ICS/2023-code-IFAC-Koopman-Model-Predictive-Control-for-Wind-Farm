function testReadInAndVisualizeKoopmanSysId

%% Set parameters
addpath(genpath(pwd))
filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat';
% filenameVal = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat';

noStates = 6; %
n = noStates;
poly = 0; % switch used for 12 states

NTMstr = '';
Vinf1 = 8; %Vinf1 = utemp(1,1,t0);
Vstr = '';
yawmode = 0;

t0 = 240;

percentTrain = 0.7;

% Get X and Y coordinates of the turbines
Wp.turbine = struct(...
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747]... % Y-coordinates of turbines (m)
    );
Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );

%% Set path to data directories
codedir = mfilename('fullpath');
parentdir = fileparts(fileparts(codedir));

dataOutDir = fullfile(parentdir,'DataInOutWfSim');
if poly == 1 && yawmode == 0
    dirdmdName = 'eDMDresults_UasOutput_poly';
elseif poly == 0 && yawmode == 0
    dirdmdName = 'eDMDresults_UasOutput';
elseif poly == 1
    dirdmdName = 'eDMDresults_UasOutput_MIMO_poly';
else
    dirdmdName = 'eDMDresults_UasOutput_MIMO';
end
dmddir = fullfile(dataOutDir,dirdmdName);

folderName = strrep(sprintf('Vin_Vinf%2.1f%s%s_states%02d',Vinf1,Vstr,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); % stateName
dirFig = fullfile(dmddir,folderName); % this depends on vinf

%% Load data: System identified model: load(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
load(fullfile(dirFig,fileName),'Koop');
K = Koop.K'; % sys_red0 = ss(K(1:n,1:n), K(1:n,n+1:tend),[eye(2),zeros(2,4)],zeros(2,3),1); step(sys_red0)

%% Load data: Simulation run
tmpId = loadSimulationRun(parentdir,filenameId);
[Inputs, Inputs_val, Outputs, Outputs_val, Deterministic,Deterministic_val,Input_u,Inputvalid_u]  = readSimulationRunData(tmpId,yawmode,percentTrain,t0,Wp); %#ok<ASGLU> 

%% Generate prediction data
InputsAll = [Inputs; Input_u];
[ysim,FITje] = eDMDsimulate(Deterministic,InputsAll,K,poly,n,yawmode);

InputsAll_val = [Inputs_val; Inputvalid_u];
[ysim_val,FITje_val] = eDMDsimulate(Deterministic_val,InputsAll_val,K,poly,n,yawmode);

%% Plot the data
strVal = 'Id.'; Vin = 0;
plotEDMDinputsEffWind(ysim,InputsAll,Deterministic, FITje,dirFig,Vin,n,strVal);
strVal = 'Val.';
plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,n,strVal);

end


function tmpId = loadSimulationRun(parentdir,filenameId)
% Load data: Simulation data
DataIn = fullfile(parentdir,'DataT2OLWFSim'); %,...
sepStr = strfind(filenameId,'_');
if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

% load identification data
dirName = fullfile(DataIn, [filenameId(1:sepStr-1),'_OL_Ct'], filenameId);
tmpId = load(dirName,'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', ...
    'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v');
end

function [Inputs, Inputs_val, Outputs, Outputs_val, Deterministic,Deterministic_val,Input_u,Inputvalid_u,Vinf] = readSimulationRunData(tmpId,yawmode,percentTrain,t0,Wp)
% get the control input, the outputs and deterministic states as well as the
% wind in front of the first turbine

fieldCell = fieldnames(tmpId);
for idx = 1:length(fieldCell)
    tmp = [tmpId.(fieldCell{idx})]; %#ok<NASGU>
    eval([fieldCell{idx},'= tmp;']);
end

tident = floor(length(Ur1)*percentTrain);
tval = tident+t0;

nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx)-1;
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);

kk = length(u)/Wp.mesh.Ny;

utemp = reshape(u,Wp.mesh.Nx,Wp.mesh.Ny,kk);
utempT = squeeze(utemp(nTx,nTy,:));
Input_u  = utempT(t0:tident)'; %#ok<*USENS>
Inputvalid_u = utempT(tval:end)';

VinfAll = squeeze(utemp(1,1,:));
Vinf = VinfAll(tval:end)';

tend = length(Ct1);

if yawmode == 0
    Inputs = [Ct1(t0:tident);Ct2(t0:tident)];
    Inputs_val = [Ct1(tval:tend);Ct2(tval:tend)];
else
    Inputs = [Ct1(t0:tident);Ct2(t0:tident); phi1(t0:tident)];
    Inputs_val = [Ct1(tval:tend);Ct2(tval:tend); phi1(tval:tend)];
end

Outputs = [PT1(t0:tident);PT2(t0:tident);FT1(t0:tident);FT2(t0:tident)];
Outputs_val = [PT1(tval:tend);PT2(tval:tend);FT1(tval:tend);FT2(tval:tend)];

Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
Deterministic_val = [Ur1(tval:tend); Ur2(tval:tend)];
end

function [ysim,FITje] = eDMDsimulate(Deterministic,InputsAll,K,poly,n,yawmode)
% eDMD_RTUpdate generates simulated output from the states

if nargin < 6,  yawmode = 0; end

% Size states, inputs, outputs
ny = 2; % number outputs:  Ur1,Ur2 (no lifted states)
Ky = K(1:ny,:); % Get the relevant rows of the Koopman matrix
structPrev.Ur1_prev1 = Deterministic(1,1); %psi_xk_1(1,1);
structPrev.Ur2_prev1 = Deterministic(2,1); %psi_xk_1(2,1);
structPrev.dUr1_prev1 = 0;
structPrev.dUr2_prev1 = 0;
structPrev.M1(1) = structPrev.Ur1_prev1; % Moving mean
structPrev.M2(1) = structPrev.Ur2_prev1;
structPrev.k = 1;

% Use non-linear Koopman functions based on deterministic states
xsim = nan(size(Ky,2),size(Deterministic,2));
xsim(:,1) = K_psi([Deterministic(:,1);InputsAll(:,1)],poly,n,structPrev,yawmode);
for idx = 1: length(Deterministic)-1
    xsim(1:ny,idx+1) = Ky * xsim(:,idx);
    [xsim(:,idx+1),~,structPrev] = K_psi([xsim(1:ny,idx+1);InputsAll(:,idx+1)], poly,n,structPrev,yawmode);
end
ysim = xsim(1:ny,:)';
FITje = max(diag(100*(eye(size(Deterministic,1))-cov(Deterministic'-ysim)./cov(Deterministic'))),0);
end

