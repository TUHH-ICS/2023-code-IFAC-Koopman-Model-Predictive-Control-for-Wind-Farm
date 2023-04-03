function mpc = MPCinitKoopman(sol,Wp,NN,refstairs,DataKoopman,KoopmanStates,R,PolyLiftingFunction,S)
% controller parameters
mpc.Nh         = 10;                              % prediction horizon
mpc.N          = Wp.turbine.N;                    % number of turbines
mpc.nu         = 2*Wp.turbine.N-1;                % number of inputs
mpc.Q          = 1e-4*eye(mpc.Nh);                % weigth on tracking
mpc.R          = R*eye(mpc.Nh*(Wp.turbine.N+1));  % weigth on control signal 1e9
mpc.duc        = 1e-1;                            % limitation on du/dt
mpc.um         = 0.2;                             % minimum CT'
mpc.uM         = 2;                               % maximum CT'
mpc.mu_old     = zeros(2*mpc.Nh*Wp.turbine.N,1);  % input matrix used in matsys function.
mpc.Isel       = kron(eye(mpc.Nh),ones(1,Wp.turbine.N));
mpc.Csel       = kron(eye(mpc.Nh*Wp.turbine.N),[0,1,0]);
mpc.C_tilde    = mpc.Isel*mpc.Csel;
mpc.Ulim_upper = repmat([mpc.uM;mpc.uM;30],mpc.Nh,1);%mpc.uM*ones(mpc.Nh*mpc.N,1);
mpc.Ulim_lower = repmat([mpc.um;mpc.um;0],mpc.Nh,1);%mpc.um*ones(mpc.Nh*mpc.N,1);
mpc.counter    = 0;

% boleans to extract desired signals from state vector
mpc.MF         = logical(repmat(repmat([1 0 0]',Wp.turbine.N,1),mpc.Nh,1));
mpc.Mf         = logical(repmat([1 0 0]',Wp.turbine.N,1));
mpc.MP         = logical(repmat(repmat([0 1 0]',Wp.turbine.N,1),mpc.Nh,1));
mpc.Mp         = logical(repmat([0 1 0]',Wp.turbine.N,1));
mpc.MU         = logical(repmat(repmat([0 0 1]',Wp.turbine.N,1),mpc.Nh,1));
mpc.Mu         = logical(repmat([0 0 1]',Wp.turbine.N,1));

% wind farm reference
N0                  = round(.1*NN/(2*Wp.sim.h));
mpc.Pref            = zeros(NN+mpc.Nh,1);
mpc.AGCdata         = S.AGCdata(:,2);
%mpc.Pgreedy         = (3.438550894184890e+07)*(Wp.site.u_Inf/12)^3/(9)*mpc.N ; % 9turb with CT'=2 in steady-state
mpc.Pgreedy         = 2.876e+06;%(7.5e+06)/(6)*mpc.N ; % from another paper6turb with CT'=2 in steady-state ; 2.875557014713873e+06 both turbines at CT'=2
if refstairs == 0
    % mpc.Pref(1:N0)      = 0.9*mpc.Pgreedy;
    % mpc.Pref(N0+1:end)  = 0.9*mpc.Pgreedy + 0.2*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);
    % mpc.Pref(1:N0)      = 0.7*mpc.Pgreedy;
    % mpc.Pref(N0+1:end)  = 0.7*mpc.Pgreedy + 0.2*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);% according to S. Boersma et al.
    mpc.Pref(1:N0)      = 0.8*mpc.Pgreedy;
    mpc.Pref(N0+1:end)  = 0.8*mpc.Pgreedy + 0.35*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);% according to S. Boersma et al.
else
    Nsep = round(NN/(5*Wp.sim.h));
    mpc.Pref(1:Nsep)      = 2e6;%0.9*mpc.Pgreedy;
    mpc.Pref(Nsep+1:2*Nsep)  = 2.2e6;%*mpc.Pgreedy+0.6e6;
    mpc.Pref(2*Nsep+1:3*Nsep)  = 2.4e6;
    mpc.Pref(3*Nsep+1:4*Nsep)  = 2.4e6;
    mpc.Pref(4*Nsep+1:end)  = 2e6;%*mpc.Pgreedy+1e6; % for palm turbines 1,1.6,1.3 greedy 1.75
end

% controller models
%KWFModel = load('Method-1_Koopman1_states4_stateNameUr1;Ur2;Ur1.^2;Ur2.^2.mat');
pathName = fileparts(pwd);
KWFModel = load(strcat(pathName,'\DataInOutWfSim\eDMDresults_UasOutput_MIMO\Vinf8dot0_states06_useval00_train7dot000000e-01\stateName_K06_P0.mat'));
%load('Method-1_Koopman1_states8_stateNameUr1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;dUr1.^3;dUr2.^3.mat');
mpc.A = KWFModel.sys_red.A; 
mpc.B = KWFModel.sys_red.B;
mpc.C = KWFModel.sys_red.C(1,:)+KWFModel.sys_red.C(2,:); % to have only P1 and P2 as output
mpc.D = KWFModel.sys_red.D(1,:)+KWFModel.sys_red.D(2,:); %outputs =[Ur1 Ur2 PT1 PT2 FT1 FT2]
mpc.nx = size(mpc.A,1);
mpc.Xinit = KWFModel.xo;
K.xprev = koopmanstateextensionWFSim(sol.u(1:2,1),PolyLiftingFunction,KoopmanStates);

mpc.xprev = K.xprev;
mpc.PolyLiftingFunction= PolyLiftingFunction;
mpc.error_norm = nan;
