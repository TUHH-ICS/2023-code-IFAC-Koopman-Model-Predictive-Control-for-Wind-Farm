%% RunWFSimDemoMeasured runs all combination of Koopman lifting function and 
% R weighting values, with and without measured values

% Clean up workspace
clear; clc; close all; 

%% User inputs

% Decide on whether stair test case (1) or defined reference is to be used
refstairs = 0;
controller = 2;

% Set R weighting values and Number of Koopman states in vectors 
RVec  = [0,1e-6,1e-4,1e-2,1,10,1e4,1e8];
KVec = [24];


%% Repeat vectors to test all combinations of R weights and Koopman operators
lenRVec = length(RVec);
lenKVec = length(KVec);

% For Koopman state 12: Run once with 12 polynomial states
PolyVec = zeros(size(KVec));
idxPolyOn = find(KVec == 12, 1, 'last');
PolyVec(idxPolyOn) = 1;

% Use repmat to repeat all Koopman states for all R values
KoopmanStatesVec = repmat(KVec,1,lenRVec); 
PolyLiftingVec = repmat(PolyVec,1,lenRVec); 
RValueVecTemp = repmat(RVec',1,lenKVec)'; 
RValueVec = RValueVecTemp(:)';

% Add combination of measured/not measured
measuredVec = [ones(size(RValueVec)), zeros(size(RValueVec))];
KoopmanStatesMeasVec = repmat(KoopmanStatesVec,1,2);
PolyLiftingMeasVec = repmat(PolyLiftingVec,1,2);
RValueVecMeasTemp = repmat(RValueVec,1,2); 

% Loop over all combinations
JR = nan(length(RValueVecMeasTemp),1);
JQ = nan(length(RValueVecMeasTemp),1);
for idx = 1 %length(RValueVecMeasTemp)
    R = RValueVecMeasTemp(idx);
    measured = measuredVec(idx);
    KoopmanStates = KoopmanStatesMeasVec(idx);
    PolyLiftingFunction = PolyLiftingMeasVec(idx);
    try
        % Actuator activity (JR) AA=%2.2e, Tracking error(JQ) TE=%2.2e[W]
        [~,JR(idx),JQ(idx)] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller); 
    catch me
        save(['Error',num2str(idx),'.mat'],'me');
    end
    
end

