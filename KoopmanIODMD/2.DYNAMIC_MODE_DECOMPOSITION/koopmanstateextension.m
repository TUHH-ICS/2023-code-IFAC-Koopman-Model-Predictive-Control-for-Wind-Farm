function [nonlobs,stateName] = koopmanstateextension(QQ_u, QQ_v, QQ_p,rho,Deterministic,PolyLiftingFunction,noStates)
if nargin < 5
    Deterministic = [];
end

if nargin < 6
    PolyLiftingFunction = 0;
end

if nargin < 7  
    noStates = 24;
end

n = 25;
%% Possible extension 1: Absolute velocity
absvel = (QQ_u.^2+QQ_v.^2).^ (1/2); %+QQ_w.^2;
absvelSqr = (QQ_u.^2+QQ_v.^2);
absvelCube = absvel .^ 3;

%% Possible extension 2: RANS  and 3
Mu = movmean(QQ_u,500);
Mv = movmean(QQ_v,500);
    M1 = nan(size(QQ_u));
    M2 = nan(size(QQ_v));
    M1(1) = Ur1(1);
    M2(1) = Ur2(1);
    for idx = 2: length(Ur1)
        idxStart = max(idx-n,1);
        idxEnd = idx;
        idxVec = idxStart:idxEnd;
        lVec = length(idxVec);
        
        M1(idx) = (lVec-1)/lVec *M1(idx-1) + 1/lVec *Ur1(idx-1);
        M2(idx) = (lVec-1)/lVec *M2(idx-1) + 1/lVec *Ur2(idx-1);
        
    end

    diffu = QQ_u - Mu;
    diffv = QQ_v - Mv;
    absveldiff = (diffu.^2+diffv.^2).^ (1/2);
    absveldiffSqr = (diffu.^2+diffv.^2);
    absveldiffCube = absveldiff.^ (3);
    Xaug1 = [QQ_v.*2;diffv.*QQ_u;QQ_v.*QQ_u];%;absveldiff.*QQ_u;absveldiffSqr.*QQ_u;];%absvel,absvelCube;absvelSqr;QQ_u.*QQ_v;diffu.*diffv;
    steadyu = QQ_u(:,100:204); %TODO: time spane should be taken as 100-1000 or moving average
    steadyv = QQ_v(:,100:204);absveldiffSqr.*QQ_u;
    steadyp = QQ_p(:,100:204);

% detu = QQ_u - mean(steadyu,2);
% detv = QQ_v - mean(steadyv,2);
% detp = QQ_p - mean(steadyp,2);

% Possible extension 2: RANS
%uv = - detu.*detv *rho;
%detuv = detu.*detv;
%sqruv = QQ_u.*QQ_u.*QQ_v;
%sqrudetv = QQ_u.*QQ_u.*detv;
%sqrvdetu = QQ_v.*QQ_v.*detu;
%Possible extension 3: Absolute velocity minus mean first values
%absveldet = detu.^2 + detv.^2 + detw.^2;
%absveldetuv = detu.^2 + detv.^2;
%% Possible extension 4
%absvelCube = absvel .^ (3/2);
%absvelCubedetuv = absveldetuv .^ (3/2);

%% Declare Non Linear Observables for states extension

%absveldet = detu.^2 + detv.^2;% + detw.^2;

%% Possible extension 4
%absvelCube = absveldet .^ (3/2);
%% Possible extension 5
if ~isempty(Deterministic)
    Ur1 = Deterministic(1,:);
    Ur2 = Deterministic(2,:);
    Ur1_prev = [Ur1(1,1),Ur1(1,1:end-1)];
    Ur2_prev = [Ur2(1,1),Ur2(1,1:end-1)];
    dUr1 = Ur1 - Ur1_prev;
    dUr2 = Ur2 - Ur2_prev;
    dUr1_prev = [dUr1(1,1),dUr1(1,1:end-1)];
    dUr2_prev = [dUr1(1,1),dUr1(1,1:end-1)];
    ddUr1 = dUr1 - dUr1_prev;
    ddUr2 = dUr2 - dUr2_prev;
    dUr1sqr = Ur1.^2 - Ur1_prev.^2;
    dUr2sqr = Ur2.^2 - Ur2_prev.^2;
    %ddUr1sqr = dUr1.^2 - dUr1_prev.^2;
    % ddUr2sqr = dUr2.^2 - dUr2_prev.^2;
    % Xaug2 =[Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;dUr1;dUr2;...
    % ddUr1;ddUr2;dUr1.*Ur2;dUr2.*Ur1;...
    % ddUr1.*Ur2;ddUr2.*Ur1];
    
    M1 = nan(size(Ur1));
    M2 = nan(size(Ur1));
    M1(1) = Ur1(1);
    M2(1) = Ur2(1);
    for idx = 2: length(Ur1)
        idxStart = max(idx-n,1);
        idxEnd = idx;
        idxVec = idxStart:idxEnd;
        lVec = length(idxVec);
        
        M1(idx) = (lVec-1)/lVec *M1(idx-1) + 1/lVec *Ur1(idx-1);
        M2(idx) = (lVec-1)/lVec *M2(idx-1) + 1/lVec *Ur2(idx-1);
        
    end
    
    %     M11 = movmean(Ur1,[25 0]);%movmean(Ur1,500);%(Ur2,500)
    %     M21 = movmean(Ur2,[25 0]);%(Ur2,500) %
    DUr1 = Ur1-M1;
    DUr2 = Ur2-M2;
    if PolyLiftingFunction == 1
        Xaug2 =[Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;Ur1.^4;Ur2.^4;...
            Ur1.^5;Ur2.^5;Ur1.^6;Ur2.^6];
        stateName = 'Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;Ur1.^4;Ur2.^4;Ur1.^5;Ur2.^5;Ur1.^6;Ur2.^6';
    else
        %dDUr1 =
        %Xaug2 = [Deterministic.^2;Deterministic.^3;...
        XaugAll =[Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;DUr1;DUr2;
            DUr1.^2;DUr2.^2;DUr1.^3;DUr2.^3;...
            DUr1.*Ur1;DUr2.*Ur2;DUr1.^2.*Ur1;DUr2.^2.*Ur2;...
            dUr1;dUr2;ddUr1;ddUr2;dUr1sqr;dUr2sqr;M1;M2]; 
        Xaug2 =  XaugAll(1:noStates,:);
        
         stateNameAll = ['Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;DUr1;DUr2;DUr1.^2;'...
             'DUr2.^2;DUr1.^3;DUr2.^3;DUr1Ur1;DUr2Ur2;DUr1.^2Ur1;DUr2.^2Ur2;dUr1;dUr2;',...
             'ddUr1;ddUr2;dUr1sqr;dUr2sqr;M1;M2'];
         stateCell = regexp(stateNameAll,';','split');
         strCellN = sprintf('''%s'';',stateCell{1:noStates});
         stateName = strCellN(2:end-1);
        %;dUr1;dUr2;...
        %ddUr1;ddUr2;dUr1sqr;dUr2sqr;M1;M2];%;ddUr1sqr;ddUr2sqr];%;dUr1.*Ur2;dUr2.*Ur1;...
        %ddUr1.*Ur2;ddUr2.*Ur1];%DUr1.^3.*Ur1;DUr2.^3.*Ur2];
        %diff1;diff2;Deterministic(1,:).*diff1.^3;Deterministic(2,:).*diff2.^3;diff1.^3;diff2.^3;Deterministic(1,:).*diff1;...
        %Deterministic(2,:).*diff2];%Deterministic(1,:).*Deterministic(2,:)
    end
end %Deterministic(1,:);Deterministic(2,:);;;diff1.^2;diff2.^2;Deterministic(1,:).*diff1.^2;Deterministic(2,:).*diff2.^2

%% Declare Non Linear Observables for states extension
nonlobs = Xaug2;


%% Unused code

