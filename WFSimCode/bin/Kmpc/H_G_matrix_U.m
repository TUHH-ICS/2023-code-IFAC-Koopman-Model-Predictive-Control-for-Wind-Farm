function [mpc] = H_G_matrix_U(Wp,mpc,sol,xinit)

%H_G_MATRIX Summary of this function goes here
%%%   Detailed explanation goes here
% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h2(Uk) < 0, at current estimate xl
% <=> min( x'H(xl)x + g(xl)x s.t. dh1(xl)x+h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with 
% Jk = E'Q E+(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss)

N = Wp.turbine.N;
Nh = mpc.Nh;
Isel   = kron(eye(mpc.Nh),ones(1,Wp.turbine.N));
Csel   = kron(eye(mpc.Nh*Wp.turbine.N),[0,1,0]);
C_tilde = Isel*Csel;
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

mpc.L_tilde = L_tilde;
mpc.S_tilde = S_tilde;
mpc.utemp = utemp;
end



