function mpc = sol_lemke_dU(mpc,Wp,xinit,sol)
      
% lemke_sol summary of this function goes here
% this function try to convert the QP into LCP and solve it using lemke
% The QP problem is 
% min f(x) = 0.5x'Hx+f'x
% subjected to (A)x>=b

Nh = mpc.Nh;
N = Wp.turbine.N;
Ulim_upper = mpc.Ulim_upper;
Ulim_lower = mpc.Ulim_lower;

% build Hessian and gradient matrix
[mpc] = H_G_matrix_dU(Wp,mpc,sol,xinit);
utemp = mpc.utemp;

% constructing the constraint matrix
Ac_ = [eye(Nh*N); -eye(Nh*N)];
bc0_ = [Ulim_lower; -Ulim_upper];
bc_ =  bc0_ - Ac_ * utemp;

% Turn QP into LCP and solve it
Hi = eye(size(mpc.H,1))/mpc.H;
M = Ac_*Hi*Ac_';
q = - Ac_*Hi*mpc.g'- bc_;

% initiating the initial condition
n = length(q);
mpc.mu_old = zeros(n,1);
temp1 = pinv(Ac_') *(mpc.H * utemp + mpc.g');
mpc.mu_old = temp1;

% lemke's Algorithm
[mu,err] = lemke(M,q,mpc.mu_old);

% if the lemke produces the error then the initila condition is calclutaed from the quadprog and asked to lemke algorithm 
if err > 0
    mpc.uTemp = myquadprog(mpc.H,mpc.g,-Ac_,(-bc0_ + Ac_ * mpc.utemp));
    disp('Lemke did not converge: Try again quadprog result');
    temp1 = pinv(Ac_') *(mpc.H *mpc.uTemp+ mpc.g');
    [mu,err] = lemke(M,q,temp1);
    mpc.counter = mpc.counter+1;  
end
if err >0
    error('Lemke did not converge');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = Hi*(Ac_'*mu-mpc.g'); %U _tedla
mpc.Xko = mpc.AA*xinit+mpc.BBt*U;
mpc.u = U(1:N:end);
mpc.U = U;
mpc.mu_old = mu;

