function [mpc] = H_G_matrix_dU(mpc,sol,xinit)

%H_G_MATRIX Summary of this function goes here
%%%   Detailed explanation goes here
% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h2(Uk) < 0, at current estimate xl
% <=> min( x'H(xl)x + g(xl)x s.t. dh1(xl)x+h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with 
% Jk = E'Q E+(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss)

Ct_prime = sol.turbine.CT_prime;
L_tilde = mpc.C_tilde*mpc.AA;
S_tilde = mpc.C_tilde*mpc.BBt;

% Jk =E'QE+(Uk-Uss)'R(Uk-Uss): E = Qsum*Qsel*(Lx+SUk-Pref % E : error  Assume U = CT_prime
mpc.utemp = kron(ones(mpc.Nh,1),Ct_prime); %diag(mpc.Ct_ss));
E_approx = (L_tilde*xinit + S_tilde* mpc.utemp - mpc.Pref(sol.k:sol.k+mpc.Nh-1));

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
deltaUvec = mpc.utemp - [Ct_prime; mpc.utemp(1:(mpc.Nh-1)*mpc.N)];
grad_dU = deltaUvec - [deltaUvec(mpc.N+1:end);zeros(mpc.N,1)];

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
R2 = [zeros((mpc.Nh-1)*mpc.N,mpc.N),-1 *eye((mpc.Nh-1)*mpc.N); zeros(mpc.N, mpc.Nh*mpc.N)];
RH = diag([2*ones((mpc.Nh-1)*mpc.N,1);ones(mpc.N,1)])+ R2 + R2';

% gradient :-
mpc.g = 2*(S_tilde'*mpc.Q*E_approx + mpc.R*grad_dU)'; 

% Hessian:-
H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R*RH);
mpc.H = 1/2*(H + H'); % for symmetry

mpc.L_tilde = L_tilde;
mpc.S_tilde = S_tilde;
end



