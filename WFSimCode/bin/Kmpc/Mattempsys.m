function mpc = Mattempsys(Wp,mpc,Nh)

Am  = mpc.A;
Bm  = mpc.B;
Cm  = mpc.C;
Dm  = mpc.D;
nx = size(Am,1);
ny = size(Cm,1);
nu = size(Bm,2);
L = zeros(nx*Nh,nx);
S = zeros(nx*Nh,nu*Nh);
%Ai = Am;
%Bi = Bm;
L(1:nx,1:nx) = Am;
S(1:nx,1:nu) = Bm;
Stemp1 = zeros(ny*Nh,nu*Nh);
Stemp1(1:ny,1:nu) = Cm*S(1:nx,1:nu);
Stemp1(1:ny,nu+1:nu+nu) = Dm;
for idx = 1:Nh-1
    L(idx*nx+1:(idx+1)*nx,1:nx) = Am*L((idx-1)*nx+1:idx*nx,1:nx);
    S(idx*nx+1:(idx+1)*nx,1:(idx+1)*nu) = [Am*S((idx-1)*nx+1:idx*nx,1:idx*nu),Bm];
    Stemp1(idx*ny+1:(idx+1)*ny,1:(idx+1)*nu) = [Cm*Am*S((idx-1)*nx+1:idx*nx,1:idx*nu),Cm*Bm];
    if idx < Nh-1
        Stemp1(idx*ny+1:(idx+1)*ny,(idx+1)*nu+1:(idx+1)*nu+nu) = Dm;
    end
end
Cc     = kron(eye(Nh),Cm);
mpc.AA  = Cc*L; % modified L
mpc.BB = Stemp1; % modified S
mpc.CC = Cc;
end
