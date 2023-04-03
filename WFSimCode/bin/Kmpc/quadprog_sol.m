function mpc = quadprog_sol(mpc,Wp,sol,xinit,isEst)

% quadprog_sol summary of this function goes here
% Solver solves for quadratic objective functions with linear constraints.
% quadprog finds a minimum for a problem specified by
% min f(x) = 0.5x'Hx+f'x
% subjected to (A)x<=b

if nargin < 5
    isEst = 0;
end

Nh = mpc.Nh;        % prediction trajectory length
N = Wp.turbine.N; %Wp.turbine.N;   % number turbines
% nu = N; %mpc.nu; %Number of turbine inputs = no WT

% Calculate the Hessian and Gradient
mpc = H_G_matrix_dU(mpc,sol,xinit);

% Set boundary conditions: Subjected to (Aco_)u <=bco_ <=> (Aco_)(uval + utemp) <=bco_
% Ac0_ = [- eye(Nh*nu); eye(Nh*nu)];
% bc0_ = [- mpc.Ulim_lower; mpc.Ulim_upper];
bc   = mpc.bc0_ - mpc.Ac0_ * mpc.utemp;

% Turning into QP and solve it%%%%%%%%%
% myoptions = optimset('Display','iter');
% [uval ,fval] = quadprog(mpc.H,mpc.g,Ac0_a,bca,[],[],[],[],[],myoptions);
[uval,fval]= myquadprog(mpc.H,mpc.g,mpc.Ac0_,bc);

%% Koopman estimate ct1,ct2,gamma -> PT
if isEst == 1
    nuK = 2*N - 1; %Ct1,Ct2...Ctn, gamma_1,... gamma_n-1
    mpcK = H_G_Kmatrix_dU(Wp,mpc,sol);
    
    % Rate limit: Diagonal and supra-diagonal
    offDiagDiff = [zeros((mpc.Nh-1)*nuK,nuK),-1*eye((mpc.Nh-1)*nuK)];
    diagDiff = eye((Nh-1)*nuK, Nh*nuK) + offDiagDiff;
    
    AcK0 = [-eye(Nh*nuK); eye(Nh*nuK)];
    AcK0d = [AcK0;- diagDiff; diagDiff];
    bc0K_ = [- repmat([mpc.Ulim_lower(1:2); 0], Nh,1); ...
        repmat([mpc.Ulim_upper(1:2); 22], Nh,1)];
    bcK = bc0K_ - AcK0 * mpcK.utemp;
    bcKd = [bcK; [- repmat([-2;-2;-1], Nh-1,1); repmat([2;2;1], Nh-1,1)]];
    
    % myoptions = optimset('Display','iter');
    % [uvalK ,fvalK] = quadprog(mpcK.H,mpcK.g,AcK0, bcK,[],[],[],[],[],myoptions);
    % [uvalK ,fvalK] = quadprog(mpcK.H,mpcK.g,AcK0d, bcKd,[],[],[],[],[],myoptions);
    [uvalK ,fvalK] = myquadprog(mpcK.H,mpcK.g,AcK0d, bcKd);
    mpc.xK = uvalK +  mpcK.utemp;
    mpc.uvalK = uvalK;
    mpc.fvalK = fvalK;
    mpc.xinitK = mpcK.xinitK;
end

mpc.x = uval + mpc.utemp;
mpc.uval = uval;
mpc.fval = fval;



