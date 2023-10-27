function [p,stateName,structPrev]  = K_psi(in1,PolyLiftingFunction,noStates,structPrev,yawmode)
%[Xaug2,stateName,structPrev] = K_psi(Deterministic,PolyLiftingFunction,noStates,structPrev)
% calculates the augmented Koopman states: [x,f(x),u]

%% Non-linear lifting function

% Different number of inputs
if nargin == 1
    PolyLiftingFunction = 0;
    noStates = 6;
    Ur1_prev1 = 0;
    Ur2_prev1 = 0;
    dUr1_prev1 = 0;
    dUr2_prev1 = 0;
    M1 = 0;
    M2 = 0;
    k = 0;
else
    % Moving mean
    Ur1_prev1 = structPrev.Ur1_prev1;
    Ur2_prev1 = structPrev.Ur2_prev1 ;
    dUr1_prev1 = structPrev.dUr1_prev1;
    dUr2_prev1 = structPrev.dUr2_prev1 ;
    M1 = structPrev.M1; % Moving mean
    M2 = structPrev.M2;
    k = structPrev.k;
end

if nargin <= 4
    yawmode = 0;
end

% Assign values in in1-vector to variables
Ur1 = in1(1,:);
Ur2 = in1(2,:);

if size(in1,1) > 2 % 3rd and 4th input thrust coefficients ct1 and ct2
    ct1 = in1(3,:);
    ct2 = in1(4,:);
end

if size(in1,1) > 5 %5th input yaw ctrl phi, 6th U1
    phi1 = in1(5,:);
    U1 = in1(6,:);
elseif size(in1,1) > 4 && yawmode %5th input yaw ctrl phi
    phi1 = in1(5,:);
elseif size(in1,1) > 4 %5th input meas. wind
    U1 = in1(5,:);
else % do nothing
end

%% Create lifting functions

Ur1_prev = [Ur1_prev1,Ur1(1,1:end-1)]; % previous states
Ur2_prev = [Ur2_prev1,Ur2(1,1:end-1)];

dUr1 = Ur1 - Ur1_prev; % Difference to previous state
dUr2 = Ur2 - Ur2_prev;

dUr1_prev = [dUr1_prev1,dUr1(1,1:end-1)];
dUr2_prev = [dUr2_prev1,dUr2(1,1:end-1)];

ddUr1 = dUr1 - dUr1_prev;
ddUr2 = dUr2 - dUr2_prev;

dUr1sqr = Ur1.^2 - Ur1_prev.^2;
dUr2sqr = Ur2.^2 - Ur2_prev.^2;

n = 25;
lVec = min(n,k);
M1 = (lVec-1)/lVec *M1 + 1/lVec *Ur1_prev;
M2 = (lVec-1)/lVec *M2 + 1/lVec *Ur2_prev;

DUr1 = Ur1 - M1;
DUr2 = Ur2 - M2;

if PolyLiftingFunction
    Xaug2 =[Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;Ur1.^4;Ur2.^4;...
        Ur1.^5;Ur2.^5;Ur1.^6;Ur2.^6];
    stateName = 'Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;Ur1.^4;Ur2.^4;Ur1.^5;Ur2.^5;Ur1.^6;Ur2.^6';
else
    XaugAll =[ Ur1; Ur2; Ur1.^2;Ur2.^2; Ur1.^3; Ur2.^3; %state 1 to 6
        DUr1; DUr2; DUr1.^2; DUr2.^2; M1; M2; %state 7 to 12
        DUr1.*Ur1; DUr2.*Ur2; % state 13 to 14 : Multiplications
        DUr1.^3; DUr2.^3;...% state 15 to 16: Cubes
        DUr1.^2.*Ur1; DUr2.^2.*Ur2;....% state 17 to 18: Squares times Cubes
        dUr1;dUr2;ddUr1;ddUr2; dUr1sqr; dUr2sqr];
    Xaug2 =  XaugAll(1:noStates,:);
    
    stateNameAll = ['Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;'... %state 1 to 6
        'DUr1;DUr2;DUr1.^2;DUr2.^2;M1;M2;'... % %state 7 to 12 diff. to mean and mean
        'DUr1Ur1;DUr2Ur2;'...% state 13 to 14 : Multiplications
        'DUr1.^3;DUr2.^3;'... % state 15 to 16: Cubes
        'DUr1.^2Ur1;DUr2.^2Ur2;',...
        'dUr1;dUr2;ddUr1;ddUr2;dUr1sqr;dUr2sqr;'];
    stateCell = regexp(stateNameAll,';','split');
    strCellN = sprintf(' %s;',stateCell{1:noStates});
    stateName = strCellN(2:end-1);
end

% Update struct
structPrev.Ur1_prev1 = Ur1;
structPrev.Ur2_prev1 = Ur2;
structPrev.dUr1_prev1 = dUr1;
structPrev.dUr2_prev1 = dUr2;
structPrev.M1 = M1;
structPrev.M2 = M2;
structPrev.k = k + 1;

if size(in1,1) > 5
    p = [Xaug2; ct1; ct2; phi1; U1];
    stateName = [stateName,';ct1;ct2;phi1;U1'];
elseif size(in1,1) > 4 && yawmode
    p = [Xaug2; ct1; ct2; phi1];
    stateName = [stateName,';ct1;ct2;phi1'];
elseif size(in1,1) > 4
    p = [Xaug2; ct1; ct2; U1];
    stateName = [stateName,';ct1;ct2;U1'];
elseif size(in1,1) > 2
    p = [Xaug2; ct1; ct2];
    stateName = [stateName,';ct1;ct2'];
else
    p = Xaug2;
end

%% unused states
% M11 = movmean(Ur1,[25 0]);%movmean(Ur1,500);%(Ur2,500)
% M21 = movmean(Ur2,[25 0]);%(Ur2,500) %

%;dUr1;dUr2;...
%ddUr1;ddUr2;dUr1sqr;dUr2sqr;M1;M2];%;ddUr1sqr;ddUr2sqr];%;dUr1.*Ur2;dUr2.*Ur1;...
%ddUr1.*Ur2;ddUr2.*Ur1];%DUr1.^3.*Ur1;DUr2.^3.*Ur2];
%diff1;diff2;Deterministic(1,:).*diff1.^3;Deterministic(2,:).*diff2.^3;diff1.^3;diff2.^3;Deterministic(1,:).*diff1;...
%Deterministic(2,:).*diff2];%Deterministic(1,:).*Deterministic(2,:)



