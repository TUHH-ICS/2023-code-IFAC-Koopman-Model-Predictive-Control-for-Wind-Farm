function [CT_prime,phi,mpc] = qLPVMPCcontroller_test1_dU(sol,Wp,mpc)

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
        
        % cons = [];
        
        % define decision variables for the windfarm
        U    = sdpvar(Wp.turbine.N*mpc.Nh,1);
        
        % build wind farm model
        mpc  = wfmodel(sol,Wp,mpc,nlIdx,mpc.Nh);
        
        % build matrices horizon
        mpc  = matsys_v(Wp,mpc,mpc.Nh,sol,xinit);
        
        % prepare
        X    = mpc.AA*xinit + mpc.BBt*U ;
        Y    = mpc.CC*X;
        
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
        
        E    = mpc.Pref(sol.k:sol.k+mpc.Nh-1)-sum(P)';
        
        % cons = [cons, mpc.um <= U <= mpc.uM]; %#ok<AGROW> little extra time
        cons = mpc.um <= U <= mpc.uM; % cons was initialized as empty above
        
        dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];%CT
        
        cost = E'*mpc.Q*E + dU'*mpc.R*dU;
        
        ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
        
        optimize(cons,cost,ops);
        
        Uopt(:,nlIdx) = value(U); %value(Y(mpc.MU))
        Popt  = value(P);
        
        temp  = (repmat(Popt(:,nlIdx),mpc.Nh,1)./(mpc.cp(1).*Uopt(:,nlIdx))).^(1/3);
        
        mpc.V = reshape(temp,Wp.turbine.N,mpc.Nh); % rotor-averaged wind speed in the horizons
        costsolver = value(cost);
        
        %%
        Isel   = kron(eye(mpc.Nh),ones(1,Wp.turbine.N));
        Csel   = kron(eye(mpc.Nh*Wp.turbine.N),[0 1 0]);
        C_tilde = Isel*Csel;
        L_tilde = C_tilde*mpc.AA;
        S_tilde = C_tilde*mpc.BBt;
        
        % cost function Jk =E'QE+(dUk)'R(dUk-Uss) with  E = Qsum*Qsel*(Lx+SUk-Pref)
        utemp = kron(ones(mpc.Nh,1),sol.turbine.CT_prime); %diag(mpc.Ct_ss));
        E_approx = (L_tilde*xinit + S_tilde* utemp - mpc.Pref(sol.k:sol.k+mpc.Nh-1));
        
        % calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
        deltaUvec = utemp - [sol.turbine.CT_prime; utemp(1:(mpc.Nh-1)*Wp.turbine.N)];
        grad_dU = deltaUvec - [deltaUvec(2:Nh*N);0];
        
        % calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
        R2 = [zeros(Nh*N-1,1),-1 *eye(Nh*N-1); zeros(1, Nh*N)];
        RH = diag([2*ones(Nh*N-1,1);1])+ R2 + R2';
        
        % gradient :-
        mpc.g = 2*(S_tilde'*mpc.Q*E_approx + mpc.R*grad_dU);
        % Hessian:-
        H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R*RH);
        mpc.H = 1/2*(H + H'); % for symmetry
        
        Nh =mpc.Nh;
        utemp = mpc.utemp;
        N = Wp.turbine.N;% number turbines
        Ulim_upper = mpc.uM*ones(Nh*N,1);
        Ulim_lower = mpc.um*ones(Nh*N,1);
        
        % %%%%%% Turn QP into LCP and solve it%%%%%%%%%
        Ac_ = -[eye(Nh*N); -eye(Nh*N)];
        bc_ = -[Ulim_lower; - Ulim_upper] - Ac_ * utemp;
        [uval ,fval] = quadprog(mpc.H,mpc.g,Ac_,bc_);
        aValue_approx = uval + utemp;
        
        X    = mpc.AA*xinit + mpc.BBt*utemp ;
        Y    = mpc.CC*X;
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
        Etemp    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        
        d_utemp   = [utemp(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            utemp(Wp.turbine.N+1:end) - utemp(1:end-Wp.turbine.N)];%CT
        
        % f(x1) = f(x0) + g(x0)*(x1-x0) + 0.5 H(x0) *(x1-x0)^2 = f(x0) +fval
        costValue02 = Etemp'*mpc.Q*Etemp + d_utemp'*mpc.R*d_utemp + fval;
        
        % f(aValue_approx)
        X    = mpc.AA*xinit + mpc.BBt*aValue_approx;
        Y    = mpc.CC*X;
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
        Etemp1    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];%CT
        
        costValue03 = Etemp1'*mpc.Q*Etemp1 + d_utemp2'*mpc.R*d_utemp2;
        
        %% using Lemke program
        
        
        mpc = lemke_sol(mpc,Wp,xinit);
        
        aValue_approxLemke = mpc.U + mpc.utemp;%Hi*(Ac_'*mu - mpc.g')
        
        X    = mpc.AA*xinit + mpc.BBt*aValue_approxLemke;
        Y    = mpc.CC*X;%yopt
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        d_utemp2   = [aValue_approxLemke(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approxLemke(Wp.turbine.N+1:end) - aValue_approxLemke(1:end-Wp.turbine.N)];
        costValue_Lemke = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;
        %costValue_Lemke = Etemp_L'*mpc.Q*Etemp_L + aValue_approxLemke'*mpc.R*aValue_approxLemke;
        
        
    end
    
    %% Assign the decision variables
    Yopt          = value(Y);
    Uopt          = Yopt(mpc.MU);
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    
    tempplot =value(U);
    figure;
    for idxPl = 1:9,subplot(3,3,idxPl); plot(tempplot(idxPl:9:end)),...
            hold on;
        plot(aValue_approxLemke(idxPl:9:end),'--'); axis tight; grid on;...
            if idxPl == 1
            title(sprintf('Time step %d', sol.k))
            end
    end
    
end
