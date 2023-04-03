%% Toy example
clear; close all; clc;
addpath(genpath(pwd))

N = 9; Wp.turbine.N = N; % Number windturbines
Nh = 10; mpc.Nh = Nh; % Number time steps prediction trajectory
noRun = 10;

mpc.A = 0.1*eye(N); mpc.B = 0.9*eye(N); mpc.C = 1*eye(N); xinit = 0*ones(N,1);

mpc.MP = true(mpc.Nh,length(xinit));
sol.k = 1;

mpc.Pref = ones(100,1);
mpc.um = 0.2;
mpc.uM = 2;
mpc.Q = eye(mpc.Nh);
mpc.R = eye(mpc.Nh*N);

sol.turbine.CT_prime = mpc.um*ones(N,1)*4; %mpc.uM;

%% Calculate input U with optimize
% Build matrices horizon, define decision variables. prepare states,
% outputs, power output, tracking error

% build matrices horizon
mpc  = matsys(Wp,mpc);
U = sdpvar(Wp.turbine.N*mpc.Nh,1); %define decision variables for windfarm

% prepare states, outputs, power output, tracking error
X    = mpc.AA*xinit + mpc.BBt*U ;
Y    = mpc.CC*X;
P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
E    = sum(P,1)' - mpc.Pref(sol.k:sol.k+mpc.Nh-1);
cons = mpc.um <= U <= mpc.uM; % set constraints

% prepare vector with change in input and cost
dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];%CT
cost = E'*mpc.Q*E + dU'*mpc.R*dU;

ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);

% run optimize
tic;
for idx = 1: noRun
    optimize(cons,cost,ops);
end
toc_opt = toc;

anU = value(U);
E_num = value(E);

costValue = value(cost);
d_anU   = [ anU(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
    anU(Wp.turbine.N+1:end)-anU(1:end-Wp.turbine.N)];%CT

costValue01 = value(E)'*mpc.Q*value(E)+ d_anU'*mpc.R*d_anU;

% solution with optimize function
Xnum    = mpc.AA*xinit + mpc.BBt*anU ;
Ynum    = mpc.CC*Xnum;

%% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h1(Uk) < 0, at current estimate xl
% <=> min( xtilde'H(xl)xtilde + g(xl)xtilde
% s.t. dh1(xl)xtilde +h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with
% Jk = E'Q E+(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss)

Isel   = kron(eye(mpc.Nh),ones(1,Wp.turbine.N));
Csel   = kron(eye(mpc.Nh*Wp.turbine.N),1);
C_tilde = Isel*Csel;

tic;
for idx = 1: noRun
    
    L_tilde = C_tilde*mpc.AA;
    S_tilde = C_tilde*mpc.BBt;
    
    % cost function Jk =E'QE+(dUk)'R(dUk-Uss) with  E = Qsum*Qsel*(Lx+SUk-Pref)
    utemp = kron(ones(mpc.Nh,1),sol.turbine.CT_prime); %diag(mpc.Ct_ss));
    E_approx = (L_tilde*xinit + S_tilde* utemp - mpc.Pref(sol.k:sol.k+mpc.Nh-1));
    
    % calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
    deltaUvec = utemp - [sol.turbine.CT_prime; utemp(1:(N*(mpc.Nh-1)))];
    grad_dU = deltaUvec - [deltaUvec(N+1:end);zeros(N,1)];
    
    % calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
    R2 = [zeros(N*(Nh-1),N),- eye(N*(Nh-1)); zeros(N,N*Nh)];
    RH = diag([2*ones(N*(Nh-1),1);ones(N,1)])+ R2 + R2';
    
    % gradient :-
    mpc.g = 2*(S_tilde'*mpc.Q*E_approx + mpc.R*grad_dU);
    % Hessian:-
    H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R*RH);
    mpc.H = 1/2*(H + H'); % for symmetry
    
    % Limits
    Ulim_upper = mpc.uM*ones(Nh*N,1);
    Ulim_lower = mpc.um*ones(Nh*N,1);
    
    Ac_ = -[eye(Nh*N); -eye(Nh*N)];
    bc_ = -[Ulim_lower; - Ulim_upper] - Ac_ * utemp;
end
toc_prepQuadprog = toc;

%% Quadprog
tic;
for idx = 1: noRun
    [uval ,fval] = myquadprog(mpc.H,mpc.g,Ac_,bc_);
end
toc_quadprog = toc;
aValue_approx = uval + utemp;

X    = mpc.AA*xinit + mpc.BBt*utemp;
Y    = mpc.CC*X;
P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
Etemp = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P,1)';

d_utemp   = [utemp(1:Wp.turbine.N) - sol.turbine.CT_prime ; ...
    utemp(Wp.turbine.N+1:end) - utemp(1:end-Wp.turbine.N)];%CT

% f(x1) = f(x0) + g(x0)*(x1-x0) + 0.5 H(x0) *(x1-x0)^2 = f(x0) +fval
costValue02 = Etemp'*mpc.Q*Etemp + d_utemp'*mpc.R*d_utemp + fval;

% f(aValue_approx)
X    = mpc.AA*xinit + mpc.BBt*aValue_approx;
Y    = mpc.CC*X;
P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
Etemp1    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P,1)';
d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
    aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];%CT

costValue03 = Etemp1'*mpc.Q*Etemp1 + d_utemp2'*mpc.R*d_utemp2;

% LCP and Lemke's algorithm
Ac_ = [eye(Nh*N); -eye(Nh*N)];
bc_ = [Ulim_lower; - Ulim_upper] - Ac_ * utemp;

tic;
for idx = 1: noRun
    Hi = eye(size(mpc.H,1))/mpc.H;
    M = Ac_*Hi*Ac_';
    q = - Ac_*Hi*mpc.g - bc_;
    
    % n = length(q);
    % mpc.mu_old = zeros(n,1);
    
    temp1 = pinv(Ac_') *(mpc.H * utemp + mpc.g);
    mpc.mu_old = temp1;
    
    [mu,err] = lemke(M,q,mpc.mu_old);
end

toc_lemke = toc;
if err >0
    error('Lemke did not converge');
end

aValue_approxLemke = Hi*(Ac_'*mu - mpc.g) + utemp;

X    = mpc.AA*xinit + mpc.BBt*aValue_approxLemke;
Y    = mpc.CC*X;
P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;%Power output
Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P,1)';

d_utemp3   = [aValue_approxLemke(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
    aValue_approxLemke(Wp.turbine.N+1:end) - aValue_approxLemke(1:end-Wp.turbine.N)];%CT

costValue_Lemke = Etemp_L'*mpc.Q*Etemp_L + d_utemp3'*mpc.R*d_utemp3;

%% Diagnostic
fprintf('Cost optimize %2.3f, quadprog %2.3f, Lemke''s %2.3f\n',...
costValue, costValue03, costValue_Lemke);
fprintf('Runtime %d runs: optimize %2.3f, quadprog %2.3f, Lemke''s %2.3f\n',...
noRun, toc_opt, toc_quadprog + toc_prepQuadprog, toc_lemke + toc_prepQuadprog);
isSqrNo = mod(sqrt(N),1) == 0;



