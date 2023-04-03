function [CT_prime,phi,mpc] = qLPVMPCcontroller_U(sol,Wp,mpc)

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
        
        % build wind farm model
        mpc  = wfmodel(sol,Wp,mpc,nlIdx,mpc.Nh);
        
        % build matrices horizon
        mpc = Matsys(Wp,mpc,mpc.Nh);
        
        
        
        mpc = sol_lemke_U(mpc,Wp,xinit,sol);
        
        approx_value = mpc.U + mpc.utemp; %Hi*(Ac_'*mu - mpc.g')

        X    = mpc.AA*xinit + mpc.BBt*approx_value;
        Y    = mpc.CC*X;%yopt
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';

        costValue_Lemke = Etemp_L'*mpc.Q*Etemp_L + approx_value'*mpc.R*approx_value;

   
    end
    
    %% Assign the decision variables
    Yopt          = value(Y);
    Uopt          = Yopt(mpc.MU);
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    %tempplot =value(U);
%     figure; 
%     for idxPl = 1:9,subplot(3,3,idxPl); %plot(tempplot(idxPl:9:end)),...
%         %hold on; 
%         plot(approx_value(idxPl:9:end),'--'); axis tight; grid on;...
%         if idxPl == 1 
%             title(sprintf('Time step %d', sol.k)) 
%         end
%     end

end
