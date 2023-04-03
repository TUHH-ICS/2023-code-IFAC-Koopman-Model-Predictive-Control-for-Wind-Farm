function [sys_red,xo]  = eDMD_RTUpdate(states,Inputs,dt,stateName)
%[sys_red,FITje,U,S,V,X,X_p,x]
% eDMD_RTUpdate builds a state space model from the states
% and input gathered in the simulation 

X      = states(:,1:end-1); % States
X_p    = states(:,2:end); % States delayed
inp    = Inputs(:,1:end-1);

nx = size(X,1); % number states: Ur1,Ur2 plus lifted states
nu = size(inp,1); % number inputs:  Ur1,Ur2 plus lifted states
ny = 2; % number inputs:  Ur1,Ur2 plus lifted states

all1 = X_p * pinv([X;inp]);% Matrix Aall: X_p = [A,B;C,D] * [X;inp]
approxA = all1(1:nx,1:nx);
approxB = all1(1:nx,nx+1:end);
approxC = [eye(ny),zeros(ny,nx-ny)];
approxD = zeros(ny,nu);

sys_red = ss(approxA,approxB,approxC,approxD,dt);

sys_red.StateName = strsplit(stateName,';');
sys_red.OutputName = {'Ur1';'Ur2'};
xo = X(:,1);%for States Ur1,Ur2 and lifted states

