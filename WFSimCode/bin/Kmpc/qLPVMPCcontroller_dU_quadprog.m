function [CT_prime,phi,mpc] = qLPVMPCcontroller_dU_quadprog(sol,Wp,mpc,K,measured)

% for now: persistent vector (we assume that the effective wind speed of
% the last time steps is known (from power, cp and control settings)
persistent vecLastUr
persistent vecLastUrEst
persistent vecCT_prime
persistent P1

persistent v;

lenVecUr = 400;
nHinf = 0; %min(10, (sol.k-2));

if isempty(vecLastUr)
    vecLastUr = NaN(2,lenVecUr);
    vecLastUrEst = NaN(2,lenVecUr);
    vecCT_prime = NaN(2,lenVecUr);
    P1 = 0;
    v = 0.2;
end


% initialize controller
mpc = MPCinit_ctrl(sol,Wp,mpc);

% check if the Koopman model should be updated
enoughData = ~any(isnan(vecLastUr(:))); % alternatively, check for sol.k > length(vecLastUr)

% Check if inf-norm is larger then v data
error_norm = norm(vecLastUr(:,end-nHinf:end)  - vecLastUrEst(:,end-nHinf:end),'inf');

% if norm(PHX([1 2 3 5 6 7])-KPHX([1 2 3 5 6 7]),'inf') > v
% error_norm = norm(Measured_X_data(error_norm_channels)-Predicted_X_from_data(error_norm_channels),'inf');
% add = [error_norm; 0];    % Output if data was added to the library
if enoughData && error_norm > v && 0
    P1 = P1 + 1; % Number of library updates
    [states, stateName] = koopmanstateextensionWFSim(vecLastUr,mpc.PolyLiftingFunction,size(K.A,1));
    
    X      = states(:,1:end-1); % States
    X_p    = states(:,2:end); % States delayed
    inp    = vecCT_prime(:,1:end-1);
    
    nx = size(X,1); % number states: Ur1,Ur2 plus lifted states
    nu = size(inp,1); % number inputs:  Ur1,Ur2 plus lifted states
    ny = 2; % number inputs:  Ur1,Ur2 plus lifted states
    approxC = [eye(ny),zeros(ny,nx-ny)];
    approxD = zeros(ny,nu);
    
    %all0 = X_p * pinv([X;inp]);% Matrix Aall: X_p = [A,B] * [X;inp]
    Psi = [X;inp];
    K.GKoop = 1/P1*((P1-1)* K.GKoop + Psi* Psi');% Update G
    K.AKoop = 1/P1*((P1-1)* K.AKoop + X* Psi'); % Update  A
    
    all1 = K.AKoop*pinv(K.GKoop); %X_p * pinv([X;inp]);%for States Ur1,Ur2 and lifted states
    
    approxA = all1(1:nx,1:nx);
    approxB = all1(1:nx,nx+1:end);
    
    % sys_red1 = ss(K.A,K.B,approxC,approxD,1);
    sys_red = ss(approxA,approxB,approxC,approxD,1);
    % norm(sys_red-sys_red1,'inf')
    
    %figure; bode(sys_red); hold on; bode(sys_red1)
    
    sys_red.StateName = strsplit(stateName,';');
    sys_red.OutputName = {'Ur1';'Ur2'};
    %xo = X(:,1);%for States Ur1,Ur2 and lifted states
    
    
    %     psi_xk_1 = K_psi(x_(end-n+1:end));   % predict xk-1
    %     psi_xk = K_psi(x(end-n+1:end));      % predict xk
    %     G = 1/P*((P-1)*G + psi_xk_1*psi_xk_1'); % Update G
    %     A = 1/P*((P-1)*A + psi_xk_1*psi_xk');   % Update A
    %     K = pinv(G)*A;  % Update Koopman operator
    %     add = [error_norm;1];   % Store the error_norm and that the Koopman operator was updated
    %[K_psi, stateName] = koopmanstateextensionWFSim(vecLastUr,mpc.PolyLiftingFunction,size(K.A,1));
    %ys_red = eDMD_RTUpdate(K_psi,vecCT_prime,1,stateName);
    
    
    K.A = sys_red.A;
    K.B = sys_red.B;
    
    %v =  error_norm; %1.2;
else
    P1 = 0;
end



% solve mpc
if sol.k>=1
    
    % xinit = [Fk Pk CT_prime] in paper.
    xinit         = zeros(mpc.nx*Wp.turbine.N,1);
    xinit(mpc.Mf) = sol.turbine.force;
    xinit(mpc.Mp) = sol.turbine.power;
    xinit(mpc.Mu) = sol.turbine.CT_prime;
    
    % nl-1 is number of times the rotor-averaged wind speeds in the horizon
    % will be updated during one sample. If nl=1, the rotor-averaged wind speeds
    % are taken constant in the horizon
    nl = 2;
    % Uopt = NaN(mpc.Nh*Wp.turbine.N,nl);
    
    % Subject to (Aco_)x <= bco_ <=>
    mpc.Ac0_ = [- eye(mpc.Nh*mpc.N); eye(mpc.Nh*mpc.N)];
    mpc.bc0_ = [- mpc.Ulim_lower; mpc.Ulim_upper];
    %bc   = bc0_ - Ac0_ * mpc.utemp;
    
    for nlIdx = 1:nl
        
        % build wind farm model
        [K,mpc]  = wfmodelNL(sol,Wp,mpc,2,K);% wfmodel(sol,Wp,mpc,nlIdx,K);% nlIdx = 2 to estimate wind for Nh
        
        % build matrices horizon for measured and estimated eff. wind speed
        mpc  = Matsys(Wp,mpc,mpc.Nh);
        
        % builf matrices for power estimation
        mpc.Kp = MatsysP(K,mpc.Nh);
        
        % solve the optimization problem and calculate time
        tic
        mpcest = mpc;
        mpcest.B = mpcest.Best;
        mpcest.BBt = mpcest.BBtest;
        
        isEst = 1;
        mpcest = quadprog_sol(mpcest,Wp,sol,xinit,isEst);
        mpc = quadprog_sol(mpc,Wp,sol,xinit);
                
        mpc.displaytime = toc;
        
        %% True wind state space model of the wind farm
        aValue_approx = mpc.x; %mpc.uval + mpc.utemp;
        X    = mpc.AA*xinit + mpc.BBt*aValue_approx;
        Y    = mpc.CC*X;                                %yopt
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh);  %Power output
        
        % tracking error value
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        
        % change in input error value
        d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];
        
        % Calculate the cost function
        mpc.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;
        
        
        %% Estimated wind state space model of the wind farm with estims
        uValue_approx_est = mpcest.x; %mpcest.uval + mpcest.utemp;
        X_est    = mpcest.AA*xinit + mpcest.BBt*uValue_approx_est;
        Y_est    = mpcest.CC*X_est;                                %yopt
        P_est    = reshape(Y(mpcest.MP),Wp.turbine.N,mpcest.Nh);  %Power output
        
        % tracking error value
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpcest.Nh-1) - sum(P_est)';
        
        % change in input error value
        d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];
        
        % Calculate the cost function
        mpcest.costValue = Etemp_L'*mpcest.Q*Etemp_L + d_utemp2'*mpcest.R*d_utemp2;
        
        aValue_approx_K = mpcest.xK;
        P_est_K    = mpc.Kp.L * mpcest.xinitK + mpc.Kp.S * mpcest.xK;
        Etemp_L_K   = mpc.Pref(sol.k:sol.k+mpcest.Nh-1) - P_est_K;
        
        
    end
    
end
%% Assign the decision variables
Yopt          = Y;
Uopt          = Yopt(mpc.MU);
Yopt_est          = Y_est;
Uopt_est          = Yopt_est(mpcest.MU);

nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx)-1;
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);
u1 = sol.u(nTx,nTy);

tau =  mpc.A(3,3); %time constant
if mpc.OptYaw == 1 && Wp.turbine.N == 2 && (measured == 1 || mpc.KoopAIC == 1) && ...
    (mpc.Pref(sol.k) > mpc.Pgreedy)
    temp = findMaxYaw(u1,Wp);
    tempPhi = [tau*sol.turbine.Phi(1) + (1-tau)* temp;0];
else
    tempPhi =zeros(Wp.turbine.N,1);
end



if measured == 1
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]); % measured wind used
elseif mpc.KoopAIC
    temp          = reshape(Uopt_est,[Wp.turbine.N,mpc.Nh]);%estimated wind used
else
    
    temp0          = reshape(aValue_approx_K,[2*Wp.turbine.N-1,mpc.Nh]);%estimated wind used
    tempPhiK           = temp0(3,1); %zeros(Wp.turbine.N,1);
    tempPhi = [tau*sol.turbine.Phi(1) + (1-tau)* tempPhiK;0];
    %phi = [tempPhi;0];
    temp = tau*sol.turbine.CT_prime + (1-tau)* temp0(1:2,1);
    %temp = temp0(1:2,1);
end
CT_prime      = temp(1:2,1);              % first action horizon

phiRLim = 1;
aDiff = tempPhi(1) - sol.turbine.Phi(1);
phi = [sol.turbine.Phi(1) + min(max(aDiff,-phiRLim),phiRLim), tempPhi(2:end)];


%% Update the effective wind speed buffers
vecLastUr(:,1:end-1)= vecLastUr(:,2:end);
vecLastUrEst(:,1:end-1)= vecLastUrEst(:,2:end);
vecLastUr(:,end) = sol.turbine.Ur;
vecLastUrEst(:,end) =  sol.xprev(1:2);

vecCT_prime(:,1:end-1)= vecCT_prime(:,2:end);
vecCT_prime(:,end) = CT_prime;

%Save ctprime
mpc.error_norm = error_norm;
mpc.P1 = P1;

% if mod((sol.k-1)/100,1) == 0 %Figure for debugging
% codedir = mfilename('fullpath');
% maindir = fileparts(fileparts(fileparts(fileparts(codedir))));
% dirFig = fullfile(maindir,'CtFig');    %define main directory
% if ~isfolder(dirFig)
%     mkdir(dirFig)
% end
%     sN = 1;
%     figure;
%     for idxPl = 1:Wp.turbine.N, subplot(Wp.turbine.N,1,idxPl);
%         plot(Uopt(idxPl:Wp.turbine.N:end)), hold on;
%         plot(Uopt_est(idxPl:Wp.turbine.N:end),'--'); axis tight; grid on;
%         if idxPl == 1
%             title(sprintf('Time step %d: WT %d', sol.k,idxPl))
%         else
%             title(sprintf('WT %d', idxPl))
%         end
%
%         if mod(idxPl-1,sN) == 0
%             ylabel('c_T [-]');
%         end
%         if idxPl > sN*(sN-1)
%             xlabel('k [-]');
%         end
%         if idxPl == Wp.turbine.N
%             legend('Opt','LCP');
%         end
%     end
%     print(gcf,fullfile(dirFig,sprintf('CTtrajectory%03d',sol.k)), '-dpng');
% end



%     figure;
%     for idxPl = 1:9,subplot(3,3,idxPl); plot(tempplot(idxPl:9:end)),...
%                 hold on;
%                 plot(aValue_approxLemke(idxPl:9:end),'--'); axis tight; grid on;...
%                 if idxPl == 1
%                     title(sprintf('Time step %d', sol.k))
%                 end
%     end

end
