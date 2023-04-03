function mpc = matsys_v(Wp,mpc,sol,xinit)

Am  = mpc.A;
Bm  = mpc.B;
Bmt = mpc.Bt;
Cm  = mpc.C;
Nh = mpc.Nh;
nx = size(Am,1);
N = Wp.turbine.N;
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
           BBt(jj*nx+1:(jj+1)*nx , (ii-1)*N+1:ii*N) = AA*AiBt(:,:,ii);
           AA = AA*A;
        end   
    end
    mpc.BBt = BBt;

    C      =  Cm;
    Cc     = kron(eye(Nh),C);
    
    toc
else
%    code by Pablo
    % tic
    L = zeros(nx*Nh,nx);
    S = zeros(nx*Nh,N*Nh);
    Ai = Am;
    Bi = Bm;
    L(1:nx,1:nx) = Am;
    S(1:nx,1:N) = Bm;
    for idx = 1:Nh-1
        L(idx*nx+1:(idx+1)*nx,1:nx) = Ai*L((idx-1)*nx+1:idx*nx,1:nx);
        S(idx*nx+1:(idx+1)*nx,1:(idx+1)*N) = [Ai*S((idx-1)*nx+1:idx*nx,1:idx*N),Bi];
    end
    C      = Cm;
    Cc     = kron(eye(Nh),C);
    % toc
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

% Calculate S_tilde and L_tilde
Isel   = kron(eye(Nh),ones(1,N));
Csel   = kron(eye(Nh*N),[0 1 0]);
C_tilde = Isel*Csel;
L_tilde = C_tilde*L;
S_tilde = C_tilde*S;

% Calculate utemp
if ~isfield(mpc, 'U')
    X    = mpc.AA*xinit + mpc.BBt*repmat(sol.turbine.CT_prime,Nh,1);
    Y = mpc.CC*X;
    utemp = Y(mpc.MU);
else
    utemp = mpc.U;
end
% Ulim_upper = mpc.uM*ones(Nh*N,1);
% Ulim_lower = mpc.um*ones(Nh*N,1);
% utemp = Ulim_lower; %0.5*(Ulim_upper+Ulim_lower); 

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
deltaUvec = utemp - [sol.turbine.CT_prime; utemp(1:(N*(mpc.Nh-1)))];
grad_dU = deltaUvec - [deltaUvec(N+1:end);zeros(N,1)];

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
R2 = [zeros(N*(Nh-1),N),- eye(N*(Nh-1)); zeros(N,N*Nh)];
RH = diag([2*ones(N*(Nh-1),1);ones(N,1)])+ R2 + R2';
E1 = (L_tilde*xinit + S_tilde* utemp - mpc.Pref(sol.k:sol.k+Nh-1));

% gradient :-
mpc.g = 2*(S_tilde'*mpc.Q*E1 + mpc.R*grad_dU)';

% Hessian:-
H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R*RH);
mpc.H = 1/2*(H + H');

%if sum(mpc.mu_old) == 0 % ToDo: Is this really the best way for mu_old
Ac_ = [eye(Nh*Wp.turbine.N);-eye(Nh*Wp.turbine.N)];
temp1 = pinv(Ac_') *(mpc.H * utemp + mpc.g');

mpc.mu_old = temp1;
%end

% temp = mu_old
mpc.L_tilde = L_tilde;
mpc.S_tilde = S_tilde;
mpc.utemp = utemp;


end
