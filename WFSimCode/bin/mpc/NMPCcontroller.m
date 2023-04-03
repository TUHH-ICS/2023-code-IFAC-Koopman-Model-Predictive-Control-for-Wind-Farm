function [CT_prime,phi,mpc] = NMPCcontroller(sol,Wp,mpc)

% initialize controller
mpc = MPCinit_ctrl(sol,Wp,mpc);

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
    nl = 1;
    Uopt = NaN(mpc.Nh*Wp.turbine.N,nl);
    
    for nlIdx = 1:nl
        
        yalmip('clear');
              
        % define decision variables for the windfarm
        U    = sdpvar(Wp.turbine.N*mpc.Nh,1);
        
        % build wind farm model
        mpc  = wfmodel(sol,Wp,mpc,nlIdx,mpc.Nh);
        
        % build matrices horizon
        mpc  = matsys_v(Wp,mpc,sol,xinit);
                
 
        % prepare
        X    = mpc.AA*xinit + mpc.BBt*U ;
        Y    = mpc.CC*X;
        
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
        
        E    = mpc.Pref(sol.k:sol.k+mpc.Nh-1)-sum(P)';
        
        cons = mpc.um <= U <= mpc.uM; % cons was initialized as empty above
        
        dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];%CT
        
        cost = E'*mpc.Q*E + dU'*mpc.R*dU;
        
        ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
        
        optimize(cons,cost,ops);
        
        Uopt(:,nlIdx) = value(U); %value(Y(mpc.MU))
        Popt       = value(P);
        
        %temp       = (repmat(Popt(:,ll),mpc.Nh,1)./(mpc.cp(1).*Uopt(:,ll)) ).^(1/3);
        temp       = (repmat(Popt(:,nlIdx),mpc.Nh,1)./(mpc.cp(1).*Uopt(:,nlIdx))  ).^(1/3);
        
    end
    
    %% Assign the decision variables
    Yopt          = value(Y);
    Uopt          = Yopt(mpc.MU);
    mpc.U = Uopt;
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    
end