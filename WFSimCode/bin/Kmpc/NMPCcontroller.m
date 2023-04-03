function [CT_prime,phi,mpc] = NMPCcontroller(sol,Wp,NN)

% initialize controller
% controller parameters
Nh         = 10;                              % prediction horizon
Q          = 1e-4*eye(Nh);                   % weigth on tracking
R          = 1e9*eye(Nh*Wp.turbine.N);      % weigth on control signal
duc        = 1e-1;                            % limitation on du/dt
um         = .1;                              % minimum CT'
uM         = 2;                               % maximum CT'

% boleans to extract desired signals from state vector
MF         = logical(repmat(repmat([1 0 0]',Wp.turbine.N,1),Nh,1));
Mf         = logical(repmat([1 0 0]',Wp.turbine.N,1));
MP         = logical(repmat(repmat([0 1 0]',Wp.turbine.N,1),Nh,1));
Mp         = logical(repmat([0 1 0]',Wp.turbine.N,1));
MU         = logical(repmat(repmat([0 0 1]',Wp.turbine.N,1),Nh,1));
Mu         = logical(repmat([0 0 1]',Wp.turbine.N,1));

mpc = MPCinit(sol,Wp,NN);

% solve mpc
if sol.k>=1
    
    xinit     = zeros(mpc.nx*Wp.turbine.N,1);
    xinit(Mf) = sol.turbine.force;
    xinit(Mp) = sol.turbine.power;
    xinit(Mu) = sol.turbine.CT_prime;
    
    % nl-1 is number of times the rotor-averaged wind speeds in the horizon
    % will be updated during one sample. If nl=1, the rotor-averaged wind speeds
    % are taken constant in the horizon
    nl = 2;
    
    for ll=1:nl
        
        yalmip('clear');
        
        cons = [];
        
        % define decision variables for the windfarm
        U    = sdpvar(Wp.turbine.N*Nh,1);
        
        % build wind farm model
        mpc  = wfmodelNMPC(sol,Wp,mpc,ll,Nh);
        
        % build matrices horizon
        mpc  = matsys(Wp,mpc,Nh);
        
        X    = mpc.AA*xinit + mpc.BBt*U ;
        Y    = mpc.CC*X;
        
        P    = reshape(Y(MP),Wp.turbine.N,Nh) ;
        
        E    = mpc.Pref(sol.k:sol.k+Nh-1)-sum(P)';
        
        cons = [cons, um <= U <= uM];
        
        dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];%
        
        cost = E'*Q*E + dU'*R*dU;
        
        ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
        
        optimize(cons,cost,ops);
        
        Uopt(:,ll) = value(U); %value(Y(mpc.MU))
        Popt       = value(P);
        
        temp       = ( repmat(Popt(:,ll),Nh,1)./(mpc.cp(1).*Uopt(:,ll))  ).^(1/3);
        mpc.V      = reshape(temp,Wp.turbine.N,Nh); % rotor-averaged wind speed in the horizons
        
    end
    
    %% Assign the decision variables
    Yopt          = value(Y);
    Uopt          = Yopt(MU);
    temp          = reshape(Uopt,[Wp.turbine.N,Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    
end