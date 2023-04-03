function [mpc] = H_G_Kmatrix_dU(Wp,mpc,sol)
%H_G_MATRIX Summary of this function goes here
%%%   Detailed explanation goes here
% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h2(Uk) < 0, at current estimate xl
% <=> min( x'H(xl)x + g(xl)x s.t. dh1(xl)x+h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with 
% Jk = E'Q E+(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss)

Ct_prime = sol.turbine.CT_prime;
gamma = sol.turbine.Phi(1,:);
nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx)-1;
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);
u1 = sol.u(nTx,nTy);
uIn = [Ct_prime;gamma];
dIn = u1;
wIn = [uIn;dIn];
nuK = length(uIn);

L = mpc.Kp.L;
S = mpc.Kp.S;

nK = size(L,2)-1;
xinitUK = K_psi([sol.turbine.Ur; wIn],0,nK);
xinitK = [xinitUK(1:nK);dIn];

% Jk =E'QE+(Uk-Uss)'R(Uk-Uss): E = Qsum*Qsel*(Lx+SUk-Pref % E : error  Assume U = CT_prime
mpc.utemp = kron(ones(mpc.Nh,1),uIn); %diag(mpc.Ct_ss));
E_approx = (L*xinitK + S* mpc.utemp - mpc.Pref(sol.k:sol.k+mpc.Nh-1));

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
deltaUvec = mpc.utemp - [Ct_prime;gamma; mpc.utemp(1:(mpc.Nh-1)*nuK)];
grad_dU = deltaUvec - [deltaUvec(nuK+1:end);zeros(nuK,1)];

% calculate hessian of dU'dU: 2*[2 -1 0 0........;
%                               -1 2 -1 0 ........;
%                                0 -1 2 -1 0.......;
%                                0 0 0 -1 2 -1 0...;
%                                ...................
%                                0 0 .........-1 2 -1;
%                                0 0 ........... -1  1]
R2 = [zeros((mpc.Nh-1)*nuK,nuK),-1 *eye((mpc.Nh-1)*nuK); zeros(nuK, mpc.Nh*nuK)];
RH = diag([2*ones((mpc.Nh-1)*nuK,1);ones(nuK,1)])+ R2 + R2';

% Ct: 0.2-2: phi = 0:25 JR = mean(mean(abs(diff(CT_prime(:,time0:end)'))))/0.18 + mean(mean(abs(diff(Phi(:,time0:end)'))))/25;
mpc_R = mpc.R(1)*diag(repmat([1, 1,0.1], 1, length(grad_dU)/3));
% %mpc_R = mpc.R(1)* eye(length(grad_dU));
% % gradient :-
% mpc.g = 2*(S'*mpc.Q*E_approx + mpc_R*grad_dU)'; 
% % Hessian:-
% H = 2*(S'*mpc.Q*S + mpc_R*RH);

RU =  mpc.RU * 10^(-5)*mpc.R(1)*diag(repmat([0,0,mpc.Rvec(3)], 1, length(grad_dU)/3));
%mpc_R = mpc.R(1)* eye(length(grad_dU));
% gradient :-
mpc.g = 2*(S'*mpc.Q*E_approx + mpc_R*grad_dU + RU*mpc.utemp)'; 
% Hessian:-
H = 2*(S'*mpc.Q*S + mpc_R*RH + RU);

mpc.H = 1/2*(H + H'); % for symmetry
mpc.L = L;
mpc.S = S;
mpc.xinitK = xinitK;

