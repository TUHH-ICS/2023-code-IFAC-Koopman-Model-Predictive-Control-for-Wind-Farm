function [mpc,K] = MPCinitv(sol,Wp,NN,refstairs,DataKoopman,KoopmanStates,R,PolyLiftingFunction,S,mpc)
% MPCinitv initialize plant and controller parameters prior simulation

% controller parameters
mpc.Nh         = 10;                              % prediction horizon
mpc.N          = Wp.turbine.N;
mpc.Q          = 1e-4*eye(mpc.Nh);                % weigth on tracking
mpc.R          = R*eye(mpc.Nh*Wp.turbine.N);    % weigth on control signal 1e9
mpc.duc        = 1e-1;                            % limitation on du/dt
mpc.um         = 0.2;                             % minimum CT'
mpc.uM         = 2;                               % maximum CT'
mpc.mu_old     = zeros(2*mpc.Nh*Wp.turbine.N,1);  % input matrix used in matsys function.
mpc.Isel       = kron(eye(mpc.Nh),ones(1,Wp.turbine.N));
mpc.Csel       = kron(eye(mpc.Nh*Wp.turbine.N),[0,1,0]);
mpc.C_tilde    = mpc.Isel*mpc.Csel;
mpc.Ulim_upper = mpc.uM*ones(mpc.Nh*mpc.N,1);
mpc.Ulim_lower = mpc.um*ones(mpc.Nh*mpc.N,1);
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
mpc.Pgreedy         = 2.876e+06; %(7.5e+06)/(6)*mpc.N ; % from another paper6turb with CT'=2 in steady-state ; 2.875557014713873e+06 both turbines at CT'=2
if refstairs == -1
    tmp0 = mpc.Pref;
    tmp0(1:N0)      = 0.9*mpc.Pgreedy;
    tmp0(N0+1:end)  = 0.9*mpc.Pgreedy + 0.25*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);
    tmp = repmat(tmp0',2,1);
    tmpVec = tmp(:);
    mpc.Pref = tmpVec(1:NN+mpc.Nh);
    % mpc.Pref(1:N0)      = 0.7*mpc.Pgreedy;
    % mpc.Pref(N0+1:end)  = 0.7*mpc.Pgreedy + 0.2*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);% according to S. Boersma et al.
elseif refstairs == -2
    tmp0 = mpc.Pref;
    
    tempAGC = mpc.AGCdata;
    
    tmp0(1:N0)      = 0.9*mpc.Pgreedy;
    tmp0(N0+1:end)  = 0.9*mpc.Pgreedy + 0.2*mpc.Pgreedy*abs(tempAGC(1:NN+mpc.Nh-N0));
    
    
    phiRLim = mpc.Pgreedy*max(abs(diff(tempAGC)));
    
    
    for idx = 2:length(tmp0) % limit the yaw rate
        aDiff = tmp0(idx) - tmp0(idx-1);
        tmp0(idx) = tmp0(idx-1) + min(max(aDiff,-phiRLim),phiRLim);
    end
    
    mpc.Pref = tmp0;

    
elseif refstairs == 0
    mpc.Pref(1:N0)      = 0.8*mpc.Pgreedy;
    mpc.Pref(N0+1:end)  = 0.8*mpc.Pgreedy + 0.35*mpc.Pgreedy*mpc.AGCdata(1:NN+mpc.Nh-N0);% according to S. Boersma et al.
elseif refstairs == 10
    mpc.Pref = 3.1*10^6*ones(NN+mpc.Nh,1);
else  
    Nsep = round(NN/(5*Wp.sim.h));
    mpc.Pref(1:Nsep)      = 2e6;%0.9*mpc.Pgreedy;
    mpc.Pref(Nsep+1:2*Nsep)  = 2.2e6;%*mpc.Pgreedy+0.6e6;
    mpc.Pref(2*Nsep+1:3*Nsep)  = 2.4e6;
    mpc.Pref(3*Nsep+1:4*Nsep)  = 2.4e6;
    mpc.Pref(4*Nsep+1:end)  = 2e6;%*mpc.Pgreedy+1e6; % for palm turbines 1,1.6,1.3 greedy 1.75
end

% controller models
mpc.tau             = 5; % time constant filter CT'
[mpc.num,mpc.den]   = tfdata(c2d(tf(1,[mpc.tau 1]),Wp.sim.h),'v');% filter on force
for kk = 1:Wp.turbine.N
    mpc.a{kk}     = kron(eye(3),-mpc.den(2));
    mpc.c{kk}     = eye(3);
end

% Koopman model: directory based on name (Vinf, no. states and PolyLiftingFunction)
if length(Wp.site.u_Inf)== 1
    VinfStr = strrep(sprintf('Vinf%2.1f',Wp.site.u_Inf),'.','dot');
else
    meanV = mean(Wp.site.u_Inf);
    VinfStr = strrep(sprintf('Vinf%2.1f',meanV),'.','dot');
end
Kstr = num2str(KoopmanStates,'%02d');
koopmanDir = fullfile([DataKoopman],['Vin_',VinfStr,'_states',Kstr]);
%koopmanDir = fullfile([DataKoopman,'_MIMO'],['Vin_',VinfStr,'_states',Kstr]);
%koopmanDir = fullfile([DataKoopman],[VinfStr,'_states',Kstr]);
koopmanDirP = fullfile([DataKoopman,'_MIMO_PT'],['Vin_',VinfStr,'_states',Kstr]);
%koopmanDirP = fullfile([DataKoopman,'_MIMO_PT'],[VinfStr,'_states',Kstr]);

% if mpc.OptYaw
%     koopmanDir = fullfile([DataKoopman,'_MIMO'],['Vin_',VinfStr,'_states',Kstr]);
% end
matKoopman = sprintf('stateName_K%s_P%d.mat',Kstr,PolyLiftingFunction);
KWFModel = load(fullfile(koopmanDir,matKoopman));

KWFModelP = load(fullfile(koopmanDirP,matKoopman));

% Koopman model for Palm 2Turbine wind speed 8m\s
K.A = KWFModel.sys_red.A;
K.B = KWFModel.sys_red.B;
K.AKoop = KWFModel.Koop.A;
K.GKoop = KWFModel.Koop.G;

% initial states: current Vinf as effective wind speeds  K.xprev = KWFModel.xo;
K.xprev = koopmanstateextensionWFSim(sol.u(1:2,1),PolyLiftingFunction,KoopmanStates);
% 
% KP.A = KWFModelP.sys_red.A;
% KP.B = KWFModelP.sys_red.B;
% KP.C = KWFModelP.sys_red.C;
% KP.D = KWFModelP.sys_red.D;

K.KPsys = KWFModelP.sys_red; 
% K.AKoop = KWFModelP.Koop.A;
% K.GKoop = KWFModelP.Koop.G;

% initial states: current Vinf as effective wind speeds  K.xprev = KWFModel.xo;
K.xprev = koopmanstateextensionWFSim(sol.u(1:2,1),PolyLiftingFunction,KoopmanStates);

mpc.xprev = K.xprev;
mpc.nx = size(mpc.a{1},1);
mpc.PolyLiftingFunction= PolyLiftingFunction;
mpc.error_norm = nan;

