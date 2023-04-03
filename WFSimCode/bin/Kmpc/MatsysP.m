function K = MatsysP(K,Nh)
%Matsys sets up the stacked matrices L and S

% matrices and their sizes
Ax  = K.KPsys.A;
Bud  = K.KPsys.B;
Cx  = K.KPsys.C;
Dud  = K.KPsys.D;

% inputs: 3 control inputs (ct1, ct2,gamma0), 1 dist v1: If d > 0, dist. is
% modelled as a state
nu = 3; % size(Bud,2)-1;
nd = size(Bud,2) - nu;
nx = size(Ax,1) + nd;
ny = size(Cx,1);

B = [Bud(:,1:nu);zeros(nd,nu)];
D = Dud(:,1:nu);

if nd > 0
    A = [[Ax,Bud(:,nu+1:end)]; [zeros(1,nx-1),1]];
    C = [Cx,Dud(:,nu+1:4)];
else
    A = Ax;
    C = Cx;
end


% initialize and generate matrices
L = zeros(ny*Nh,nx);
S = zeros(ny*Nh,nu*Nh);

L(1:ny,1:nx) = C;
S(1:ny,1:nu) = D;

for idx = 1:Nh-1
    L(idx*ny+1:(idx+1)*ny,1:nx) = L((idx-1)*ny+1:idx*ny,1:nx)*A;
    S(idx*ny+1:(idx+1)*ny,1:(idx+1)*nu) = [C*A^(idx-1)*B,S((idx-1)*ny+1:idx*ny,1:idx*nu)];
end

% constant c matrix
%Cc     = kron(eye(Nh),C);

% output matices
K.L  = L;
K.S = S;
% K.BBtest = Sest;
% K.CC  = Cc;
end

