function [CT_prime,phi,mpc] = KMPCcontroller_dU_quadprog(sol,Wp,mpc)
% initialize controller
%mpc = MPCinit_ctrl(sol,Wp,mpc);

% solve mpc
if sol.k>=1
   
    % xinit = [Fk Pk CT_prime] in paper.
    %xinit         = zeros(mpc.nx,1);
    xinit         = mpc.Xinit; %TODO: We need to change initial condition 
                                % according to new data
%     xinit         = zeros(mpc.nx*Wp.turbine.N,1);
%     xinit(mpc.Mf) = sol.turbine.force;
%     xinit(mpc.Mp) = sol.turbine.power;
%     xinit(mpc.Mu) = sol.turbine.CT_prime;
    
    % nl-1 is number of times the rotor-averaged wind speeds in the horizon
    % will be updated during one sample. If nl=1, the rotor-averaged wind speeds
    % are taken constant in the horizon
    nl = 1;
    Uopt = NaN(mpc.Nh*Wp.turbine.N,nl); 
     % build matrices horizon
        mpc  = Mattempsys(Wp,mpc,mpc.Nh);
        tic
        mpc = quadprog_sol(mpc,Wp,sol,xinit);
        mpc.displaytime = toc;
        
        % state space model of the wind farm
        aValue_approx = mpc.uval + mpc.utemp; %Ct(k+1) = dCt(k)+Ct(k)
        Y    = mpc.AA*xinit + mpc.BB*aValue_approx; %yopt=LX+SU
        %Y    = mpc.CC*X;                                %yopt
        P    = Y;%reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh);  %Power output
        
        % tracking error value
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        
         % change in input error value
        d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];
        
        % Calculate the cost function
        mpc.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;


end
    
 %% Assign the decision variables
    Yopt          = Y;
    %cttemp=[aValue_approx(1:Wp.turbine.N)+sol.turbine.CT_prime ; ...
            %aValue_approx(Wp.turbine.N+1:end) + aValue_approx(1:end-Wp.turbine.N)];
    Uopt          = aValue_approx;%Yopt(mpc.MU);
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    
end