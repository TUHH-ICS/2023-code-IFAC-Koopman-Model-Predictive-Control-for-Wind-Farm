function [sol_array,Power,CT_prime,Phi,V,JR,time] = getInfoPlot(fileName,DataOut)
%getInfoPlot gets the info for the plot

tmp = load(fullfile(DataOut,fileName));
sol_array = struct;

% Control inputs for plot and eDM
CT_prime(1,:) = tmp.Ct1; CT_prime(2,:) = tmp.Ct2;
Phi(1,:) = tmp.phi1; Phi(2,:) = tmp.phi2;
time = 1:length(tmp.Ct1);

% Values open loop: Power, force and effective wind speed at rotor
Power(1,:) = tmp.PT1; Power(2,:) = tmp.PT2;
V(1,:)= tmp.Ur1; V(2,:) = tmp.Ur2;

JR = mean(mean(abs(diff(CT_prime'))));