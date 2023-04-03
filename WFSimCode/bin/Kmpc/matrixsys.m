function mpc = matrixsys(Wp,mpc,Nh,sol,xinit)

Am  = mpc.A;
Bm  = mpc.B;
Bmt = mpc.Bt;
Cm  = mpc.C;
nx = size(Am,1);
nu = Wp.turbine.N;
K=0;

if K == 1
    tic
    
    A = Am;
    Ai = A;
    AA = Ai;
    for ii = 2:Nh
        Ai = A*Ai;
        AA = [AA;Ai];
    end

    AiB = Bm;
    BB = kron(eye(Nh),AiB);
    for ii = 1:Nh-1
        AiB = A*AiB;
        BB = BB+kron(diag(ones(Nh-ii,1),-ii),AiB);
    end
    mpc.BB  = BB;

    AiBt = Bmt;
    q    = num2cell(AiBt,[1,2]);
    BBt  = blkdiag(q{:});
    for ii = 1:Nh-1
        AA = A;
        for jj = ii:Nh-1
           BBt(jj*nx+1:(jj+1)*nx , (ii-1)*nu+1:ii*nu) = AA*AiBt(:,:,ii);
           AA = AA*A;
        end   
    end
    mpc.BBt = BBt;

    C      =  Cm;
    Cc     = kron(eye(Nh),C);
    
    toc
else
%    code by Pablo
    tic
    L = zeros(nx*Nh,nx);
    S = zeros(nx*Nh,nu*Nh);
    Ai = Am;
    Bi = Bm;
    L(1:nx,1:nx) = Am;
    S(1:nx,1:nu) = Bm;
    for idx = 1:Nh-1
        L(idx*nx+1:(idx+1)*nx,1:nx) = Ai*L((idx-1)*nx+1:idx*nx,1:nx);
        S(idx*nx+1:(idx+1)*nx,1:(idx+1)*nu) = [Ai*S((idx-1)*nx+1:idx*nx,1:idx*nu),Bi];
    end
    C      = Cm;
    Cc     = kron(eye(Nh),C);
    toc
end

mpc.AA  = L;
mpc.BBt = S;
mpc.CC  = Cc;

%% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h2(Uk) < 0, at current estimate xl
% <=> min( x'H(xl)x + g(xl)x s.t. dh1(xl)x+h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with 
% Jk = E'Q E+(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss)
Isel   = kron(eye(mpc.Nh),ones(1,Wp.turbine.N));
Csel   = kron(eye(mpc.Nh*Wp.turbine.N),[0 1 0]);
C_tilde = Isel*Csel;
L_tilde = C_tilde*L;
S_tilde = C_tilde*S;
% Jk =E'QE+(Uk-Uss)'R(Uk-Uss)
% E = Qsum*Qsel*(Lx+SUk-Pref)
% E : error 
% Assume U = CT_prime
utemp = kron(ones(mpc.Nh,1),sol.turbine.CT_prime) ; 
E = (L_tilde*xinit+S_tilde* utemp - mpc.Pref(sol.k:sol.k+mpc.Nh-1));

% gradient :-
% mpc.g = 2*(S_tilde'*mpc.Q*E +...
%     [1;zeros(89,1)]'*mpc.R*(utemp - [sol.turbine.CT_prime; utemp(1:(mpc.Nh-1)*Wp.turbine.N)]))';
mpc.g = 2*(S_tilde'*mpc.Q*E +...
    *mpc.R*(utemp - [sol.turbine.CT_prime; utemp(1:(mpc.Nh-1)*Wp.turbine.N)]))';

% Hessian:-
mpc.H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R);

if sum(mpc.mu_old) == 0

    Ac_ = [eye(mpc.Nh*Wp.turbine.N);-eye(mpc.Nh*Wp.turbine.N)];
    temp1 = pinv(Ac_') *(mpc.H * utemp + mpc.g');
    
    mpc.mu_old = temp1;
end
% temp = mu_old
mpc.L_tilde =L_tilde;
mpc.S_tilde = S_tilde;
mpc.utemp = utemp;



end
