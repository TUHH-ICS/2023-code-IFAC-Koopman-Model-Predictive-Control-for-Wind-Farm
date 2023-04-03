function mpc = sol_lemke_U(mpc,Wp,xinit,sol)

% build Hessian and gradient matrix
[mpc] = H_G_matrix_U(Wp,mpc,sol,xinit);
        
Nh =mpc.Nh;
utemp = mpc.utemp;
N = Wp.turbine.N;
Ulim_upper = mpc.uM*ones(Nh*N,1);
Ulim_lower = mpc.um*ones(Nh*N,1);

%% %%%%%%Turn QP into LCP and solve it%%%%%%%%%
Ac_ = [eye(Nh*N); -eye(Nh*N)];
bc_ = [Ulim_lower; - Ulim_upper] - Ac_ * utemp;

Hi = eye(size(mpc.H,1))/mpc.H;
M = Ac_*Hi*Ac_';
q = - Ac_*Hi*mpc.g- bc_;

n = length(q);
mpc.mu_old = zeros(n,1);
temp1 = pinv(Ac_') *(mpc.H * utemp + mpc.g);
mpc.mu_old = temp1;

[mu,err] = lemke(M,q,mpc.mu_old);

if err >0
    error('Lemke did not converge');
end

U = Hi*(Ac_'*mu-mpc.g); %U _tedla


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpc.Xko = mpc.AA*xinit+mpc.BBt*U;
mpc.u = U(1:N:end);
mpc.U = U;
mpc.mu_old = mu;
