function [CT_prime,phi,mpc] = KMPCcontroller_test(sol,Wp,mpc)
% initialize controller
%mpc = MPCinit_ctrl(sol,Wp,mpc);

% solve mpc
if sol.k>=1
   
    xinit         = mpc.Xinit; %TODO: We need to change initial condition 
                                % according to new data
    nl = 1;
     yalmip('clear');
    Uopt = NaN(mpc.Nh*size(mpc.B,2),nl); 
    
     U    = sdpvar(mpc.Nh*size(mpc.B,2),nl);
     % build matrices horizon
        mpc  = Mattempsys(Wp,mpc,mpc.Nh);
        tic
        X    = mpc.AA*xinit + mpc.BB*U ;
       % Y    = mpc.CC*X + mpc.CC*U
        E    = mpc.Pref(sol.k:sol.k+mpc.Nh-1)-X;
        %Aco = [eye(mpc.Nh*(Wp.turbine.N+1));-eye(mpc.Nh*(Wp.turbine.N+1))];
        bco_min = repmat([mpc.um;mpc.um;0],mpc.Nh,1);
        bco_max = repmat([mpc.uM;mpc.uM;30],mpc.Nh,1);
        %cons = Aco*U-bco <=0;
        cons = (U >= bco_min) &( U <= bco_max); % cons was initialized as empty above
        
        dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ;U(Wp.turbine.N+1)-sol.turbine.Phi(1,:) ; U(Wp.turbine.N+2:end)-U(1:end-Wp.turbine.N-1)];%CT and gamma
        
        cost = E'*mpc.Q*E + dU'*mpc.R*dU;
        
        ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
        
        optimize(cons,cost,ops);
        mpc = quadprog_sol(mpc,Wp,sol,xinit);

        
        
        Uopt(:,1) = value(U); %value(Y(mpc.MU))
        %Xopt       = value(X);
         
        if mod((sol.k-1)/100,1) == 0
            sN = 1;
            figure;
            for idxPl = 1:Wp.turbine.N+1, subplot(Wp.turbine.N+1,1,idxPl);
                plot(Uopt(idxPl:mpc.nu:end)), hold on;
                plot(mpc.x(idxPl:mpc.nu:end),'--'); axis tight; grid on;
                if idxPl == 1
                    title(sprintf('Time step %d: WT %d', sol.k,idxPl))
                elseif idxPl == 2
                    title(sprintf('WT %d', idxPl))
                elseif idxPl == 3
                    title(sprintf('WT 1'))
                end
                
                if mod(idxPl-1,sN) == 0
                    ylabel('c_T [-]');
                end
                if idxPl > sN*(sN-1)
                    xlabel('k [-]');
                end
                if idxPl == Wp.turbine.N
                    legend('Opt','LCP');
                end
            end
            print(gcf,sprintf('CTtrajectory%03d',sol.k), '-dpng');
        end
     

end
    
 %% Assign the decision variables
  
    temp          = reshape(Uopt,[Wp.turbine.N+1,mpc.Nh]);
    %mpc.Xinit     = Xopt(:,1);
    CT_prime      = temp(1:2,1);              % first action horizon
    phi(1,1)      = temp(3,1);
    phi(2,1)      = zeros(1,1);
    mpc.Xinit     = mpc.A*xinit + mpc.B*[CT_prime;phi(1,1)];
end