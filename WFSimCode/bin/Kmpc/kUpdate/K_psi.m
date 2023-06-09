function [p,stateName,structPrev]  = K_psi(in1,PolyLiftingFunction,noStates,structPrev)
%[Xaug2,stateName,structPrev] = koopmanstateextensionRT(Deterministic,PolyLiftingFunction,noStates,structPrev)
% PHI
% P = PHI(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    09-Jan-2019 13:38:24

%% Non-linear lifting function
Ur1 = in1(1,:);
Ur2 = in1(2,:);

if size(in1,1) >= 4
    ct1 = in1(3,:);
    ct2 = in1(4,:);
end

if size(in1,1) == 6
    gamma =  in1(5,:);
    U1 = in1(6,:);

elseif size(in1,1) == 5
    U1 = in1(5,:);
end


if nargin == 1
    PolyLiftingFunction = 0;
    noStates = 6;
    Ur1_prev1 = 0;
    Ur2_prev1 = 0;
    dUr1_prev1 = 0;
    dUr2_prev1 = 0;
    M1 = 0;
    M2 = 0;
    k = 1;
elseif nargin == 3
    Ur1_prev1 = Ur1;
    Ur2_prev1 = Ur2;
    dUr1_prev1 = 0;
    dUr2_prev1 = 0;
    M1 = Ur1; % Moving mean
    M2 = Ur2;
    k = 1;
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

if size(in1,1) > 6
    p = [Xaug2; ct1; ct2; gamma1; U1];
    stateName = [stateName,';ct1;ct2;gamma;U1'];
    
elseif size(in1,1) > 4
    p = [Xaug2; ct1; ct2; U1];
    stateName = [stateName,';ct1;ct2; U1'];
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



