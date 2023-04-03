function [time, CT_prime, Phi, Power, Force, V, JR] = saveDataEDMD(Wp,sol_array,controller, ControlSetStr,DataOut)
% saveDataEDMD saves data into variables expected by main in Koopman SysID
% [time, CT_prime, Phi, Power, Force, Ur1, JR] = saveDataEDMD(Wp,...
% sol_array,controller, ControlSetStr,DataOut)


if ~nargin
    % load('sowfa_2turb_yaw_alm_turbl.mat');
    % load('sowfa_2turb_yaw_alm_combined.mat');
    load('sowfa2turb_alm_turb.mat'); %#ok<LOAD> This can be used to load a mat-file
end

% dataInMainDir = fullfile(parentDir,'DataOut2');
% matInDir = fullfile(dataInMainDir,dataDir,matDataDir);
% load(matInDir);

greedy = 0;
if greedy == 1
    Control = 'greedy';
elseif controller== 0
    Control = ControlSetStr;
else
    Control = num2str(controller);
end


fileName = [strrep(sprintf('Vinf%2.1f_%s',Wp.site.u_Inf,Control),'.','dot'),'.mat'];

% time for plot
time = cell2mat(arrayfun(@(x)(x.k), sol_array,'UniformOutput', false));

% Control inputs for plot and eDM
CT_prime = cell2mat(arrayfun(@(x)(x.turbine.CT_prime),sol_array,'UniformOutput', false));
Ct1 = CT_prime(1,:);  Ct2 = CT_prime(2,:);
Phi = cell2mat(arrayfun(@(x)(x.turbine.Phi),sol_array,'UniformOutput', false));

phi1 = Phi(1,:); phi2 = Phi(2,:);

% Values open loop: Power, force and effective wind speed at rotor
Power = cell2mat(arrayfun(@(x)(x.turbine.power),sol_array,'UniformOutput', false));
PT1 = Power(1,:); PT2 = Power(2,:);
Force = cell2mat(arrayfun(@(x)(x.turbine.force),sol_array,'UniformOutput', false));
FT1 = Force(1,:); FT2 = Force(1,:);
V = cell2mat(arrayfun(@(x)(x.turbine.Ur),sol_array,'UniformOutput', false));
Ur1 = V(1,:);  Ur2 = V(2,:);

% Values open loop: wind grid u,v and pressure grid 
u = cell2mat(arrayfun(@(x)(x.u), sol_array,'UniformOutput', false));
v = cell2mat(arrayfun(@(x)(x.v), sol_array,'UniformOutput', false));
p = cell2mat(arrayfun(@(x)(x.p), sol_array,'UniformOutput', false));

JR = mean(mean(abs(diff(CT_prime'))));


if greedy == 1
    PtotalGreedy = sum(PT1(1,end)+PT1(1,end));
else
    PtotalGreedy =0;
end
save(fullfile(DataOut,fileName),...
    'u','v','p','PT1','PT2', 'FT1','FT2','Ct1','Ct2','Ur1','Ur2','PtotalGreedy','phi1','phi2')

save(fullfile(DataOut,strrep(fileName,'.mat','_solArray.mat')),...
    'Wp','sol_array')


%load('Wfsim2Turb8SofaConfigqlmpcCL.mat')

