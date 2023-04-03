function [DataKoopman,DataOut,S] = WFSim_addpaths(controller, PolyLiftingFunction)

WFSimFolder = fileparts(mfilename('fullpath')); % Get WFSim directory
addpath(genpath(fullfile(WFSimFolder, 'libraries')));          % Add external libraries

%WFSim_addpaths add necessary paths for WFSim
addpath(fullfile(WFSimFolder,'layoutDefinitions')); % Folder with predefined wind farm layouts
addpath(fullfile(WFSimFolder,'controlDefinitions')); % Make use of a predefined timeseries of control inputs
addpath(fullfile(WFSimFolder,'solverDefinitions')); % Folder with model options, solver settings, etc.
addpath(fullfile(WFSimFolder,genpath('YALMIP-master')));% YALMIP folder for solver

binFolder = fullfile(WFSimFolder,'bin');

addpath(fullfile(binFolder,'analysis'));  % Add analysis files
addpath(fullfile(binFolder, 'core'));     % Add core files

mpcPath = fullfile(binFolder, 'mpc');
S = load(fullfile(mpcPath,'P_reference.mat')); % initialize controller settings

kmpcPath = fullfile(binFolder, 'Kmpc');
if controller == 1
    if contains(path, mpcPath)
        rmpath(kmpcPath);
    end   
    addpath(fullfile(binFolder, 'mpc')); % Add mpc files
elseif controller >1
    if contains(path, mpcPath)
        rmpath(mpcPath);
    end  
    addpath(fullfile(binFolder, 'Kmpc')); % Add mpc files
    addpath(fullfile(binFolder, 'Kmpc','kUpdate')); % Add mpc files
else % add nothing here
end

parentdir = fileparts(WFSimFolder);
dataKoopmanMain = fullfile(parentdir,'DataInOutWfSim');


if controller == 2 && ~isfolder(dataKoopmanMain)
    error('No Koopman model exists. Generate models first');
elseif controller ~= 2
    DataKoopman = '';
elseif PolyLiftingFunction == 1
    DataKoopman = fullfile(dataKoopmanMain,'eDMDresults_UasOutput_poly');
elseif controller == 3
    DataKoopman = fullfile(dataKoopmanMain,'eDMDresults_UasOutput_MIMO');
else
    DataKoopman = fullfile(dataKoopmanMain,'eDMDresults_UasOutput');
end

if controller == 2
addpath(DataKoopman)
end



DataOut = fullfile(parentdir,'DataT2OLWFSim');
if ~exist(DataOut,'dir')
    mkdir(DataOut);
end



% %% YALMIP (LMI solver)
% rootDir = fileparts(pwd);
% YALMIPdir = fullfile(rootDir,'YALMIP');
% SedumiDir = fullfile(rootDir,'sedumi');
%
% if isdir(YALMIPdir) %#ok<ISDIR> I like isdir!
%     addpath(genpath(YALMIPdir));
%     addpath(genpath(SedumiDir));
% end
