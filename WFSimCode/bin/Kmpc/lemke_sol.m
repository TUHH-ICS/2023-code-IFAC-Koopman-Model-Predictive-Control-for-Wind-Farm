function mpc = lemke_sol(mpc,Wp,xinit,sol)

% build Hessian and gradient matrix
[mpc] = H_G_matrix_dU(Wp,mpc,sol,xinit);

Nh = mpc.Nh; % prediction trajectory length
N = Wp.turbine.N; % number turbines
Ulim_upper = mpc.uM*ones(Nh*N,1);
Ulim_lower = mpc.um*ones(Nh*N,1);

%% %%%%%%Turn QP into LCP and solve it%%%%%%%%%
Ac0_ = [eye(Nh*N); -eye(Nh*N)];
bc0_ = [Ulim_lower; - Ulim_upper];
bc_ =  bc0_ - Ac0_ * mpc.utemp;

Hi = eye(size(mpc.H,1))/mpc.H;
M = Ac0_*Hi*Ac0_';
q = - Ac0_*Hi*mpc.g' - bc_;

n = length(q);
mpc.mu_old = zeros(n,1);
temp1 = pinv(Ac0_') *(mpc.H * mpc.utemp + mpc.g');
mpc.mu_old = temp1;

[mu,err] = lemke(M,q,mpc.mu_old);

uTemp = myquadprog(mpc.H,mpc.g,-Ac0_,(-bc0_ + Ac0_ * mpc.utemp));
mpc.x = uTemp + mpc.utemp;

if err > 0
    disp('Lemke did not converge: Try again quadprog result');
    temp1 = pinv(Ac0_') *(mpc.H * uTemp + mpc.g');
    % temp1 = pinv(Ac0_') *(mpc.H * Ulim_lower+ mpc.g');
    [mu,err] = lemke(M,q,temp1);
          
end

% if err > 0
%     disp('Lemke did not converge: Try again with upper limit');
%     %temp1 = pinv(Ac0_') *(mpc.H * uTemp + mpc.g');
%      temp1 = pinv(Ac0_') *(mpc.H * Ulim_upper+ mpc.g');
%     [mu,err] = lemke(M,q,temp1);
%           
% end

if err >0
  error('Lemke did not converge');
end

U = Hi*(Ac0_'*mu - mpc.g') + mpc.utemp;

  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpc.Xko = mpc.AA*xinit+mpc.BBt*U;

mpc.u = U(1:N:end);
mpc.U = U;
% mpc.mu_old = mu; % ToDo should this be included at some point?
