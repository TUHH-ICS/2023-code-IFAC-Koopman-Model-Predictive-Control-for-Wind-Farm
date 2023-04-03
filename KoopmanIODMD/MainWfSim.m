function MainWfSim(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)
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
%clc; clear; close all;
if ~nargin
    % directory for identification and validation data, yaw and pitch
    yawmode = 0; %0 for ct (pitch/torque) control, 1 for additional yaw control
    
    filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat'; NTMstr = '';
    filenameVal = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat'; NTMstr = '';

    
    % filenameId = 'Vinf8dot5_sowfa_2turb_yaw_alm_uniform'; %'Vinf8dot5_sowfa_2turb_yaw_alm_combined.mat'; %'Vinf8dot5_sowfa_2turb_yaw_alm_turbl_AllComb.mat';
    % filenameVal = 'Vinf8dot5_sowfa_2turb_yaw_alm_uniform'; %'Vinf8dot5_sowfa_2turb_yaw_alm_combined.mat';
    
    % User input
    noStates = 6; %number of Koopman states [,12,12,14,18,24]
    useVal = 1; % use extra data for validation
    percentTrain = .6;% 60% data for validation
    NTMstr = '';
    PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
    idxPolyOn = find(noStates == 12, 1, 'first');
    PolyVec(idxPolyOn) = 1;
    t0 = 1;   
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

DataIn = fullfile(parentdir,'DataT2OLWFSim');
if ~isfolder(DataIn)
    warning('DataIn directory does not exist')
    return
end

%% Loop for  'Vinf8dot5_sowfa_2turb_yaw_alm_combined2.mat'
sepStr = strfind(filenameId,'_');
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
        Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
        [states,stateName] = koopmanstateextensionWFSim(Deterministic,poly,n);
        
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
    folderName = strrep(sprintf('Vinf%2.1f%s%s_states%02d',Vinf,Vstr,NTMstr,n),'.','dot');%'Vinf%d_diffComb_states%d'
    fileName = sprintf('stateName_K%02d_P%d.mat',n,poly); %stateName);
    
    dirFig = fullfile(dirdmd,folderName); % this depends on vinf
    if ~exist(dirFig,'dir')
        mkdir(dirFig);
    end
    
    Deterministic = [Ur1(t0:tident); Ur2(t0:tident)];
    Deterministic_val = [Ur1(tval:end); Ur2(tval:end)];
    %     [sys_red,FITje,Xd,Xd_p,x,FITje_val,fig1,xo,Koop] = ...
    %         eDynamicmodedecomposition(states,statesvalid,Inputs,Outputs,...
    %         Deterministic,Deterministic_val,Inputs_val,Outputs_val,method,codedir,dirFig,dt,stateName,dirFig,yawmode);
    
    % Sys id with ct1,ct2 as inputs
    [sys_red,FITje,Xd,Xd_p,x,FITje_val,xo,Koop,...
        ysim,Inputs0,Deterministic0,ysim_val,Inputs_val0,Deterministic_val0]= ...
        eDMD_RT_Vin(states,stateName,statesvalid,poly,n,yawmode);
    strVal = 'Id.'; Vin = 0;
    plotEDMDinputsEffWind(ysim,Inputs0,Deterministic0, FITje,dirFig,Vin,n,strVal);
    strVal = 'Val.';
    plotEDMDinputsEffWind(ysim_val,Inputs_val0,Deterministic_val0,FITje_val,dirFig,Vin,n);
    
    
    save(fullfile(dirFig,fileName),'sys_red','FITje','FITje_val','dirName','xo','Koop');
    
    strVAF{idx} = sprintf('%s \t\t %d\t\t\t %d\t\t\t %d \t\t%d\n',num2str(n,'%02d'),...
        round(FITje(1)), round(FITje_val(1)),round(FITje(2)), round(FITje_val(2)));
    
end

fid = fopen(['VAF_',strrep(filenameId,'.mat',''),'.txt'],'w');
fprintf(fid,'NoK.\tT1(Id)\t\tT1(Val)\t\t T2(Id)\t\tT2(Val)\n');

for idx = 1: length(noStates)
    fprintf(fid,'%s',strVAF{idx});
end
fclose(fid);

%% Unused code
% Turbine and flow characteristics to be used
% rho = 1.20; %air density in [kg m^-3]
% D = 126.4; %Rotor Diameter used in simulations = 110 [m] %ToDo AD Check this in WFSim

