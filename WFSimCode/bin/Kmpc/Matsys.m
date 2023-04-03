function mpc = Matsys(Wp,mpc,Nh)
%Matsys sets up the stacked matrices L and S

% matrices and their sizes
A  = mpc.A;
B  = mpc.B;  
Best = mpc.Best;
C  = mpc.C;
nx = size(A,1);
nu = Wp.turbine.N;

% initialize and generate matrices
L = zeros(nx*Nh,nx);
S = zeros(nx*Nh,nu*Nh);
Sest = zeros(nx*Nh,nu*Nh);
L(1:nx,1:nx) = A;
S(1:nx,1:nu) = B;
Sest(1:nx,1:nu) = Best;
for idx = 1:Nh-1
    L(idx*nx+1:(idx+1)*nx,1:nx) = A*L((idx-1)*nx+1:idx*nx,1:nx);
    B = mpc.Bt(:,:,idx+1);
    Best = mpc.Btest(:,:,idx+1);
    S(idx*nx+1:(idx+1)*nx,1:(idx+1)*nu) = [A*S((idx-1)*nx+1:idx*nx,1:idx*nu),B]; 
    Sest(idx*nx+1:(idx+1)*nx,1:(idx+1)*nu) = [A*Sest((idx-1)*nx+1:idx*nx,1:idx*nu),Best]; 
end

% constant c matrix
Cc     = kron(eye(Nh),C);

% output matices
mpc.AA  = L;
mpc.BBt = S;
mpc.BBtest = Sest;
mpc.CC  = Cc;
end

