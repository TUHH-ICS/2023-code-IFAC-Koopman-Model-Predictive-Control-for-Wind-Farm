function mpc = matsys(Wp,mpc)

Am  = mpc.A;
Bm  = mpc.B;
% Bmt = mpc.Bt;
Cm  = mpc.C;
nx = size(Am,1);
nu = Wp.turbine.N;
Nh = mpc.Nh;

%    code by Pablo
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

Cc     = kron(eye(Nh),Cm);

mpc.AA  = L;
mpc.BBt = S;
mpc.CC  = Cc;


% K=0;

% if K == 1
%     tic
%     
%     A = Am;
%     Ai = A;
%     AA = Ai;
%     for ii = 2:Nh
%         Ai = A*Ai;
%         AA = [AA;Ai];
%     end
% 
%     AiB = Bm;
%     BB = kron(eye(Nh),AiB);
%     for ii = 1:Nh-1
%         AiB = A*AiB;
%         BB = BB+kron(diag(ones(Nh-ii,1),-ii),AiB);
%     end
%     mpc.BB  = BB;
% 
%     AiBt = Bmt;
%     q    = num2cell(AiBt,[1,2]);
%     BBt  = blkdiag(q{:});
%     for ii = 1:Nh-1
%         AA = A;
%         for jj = ii:Nh-1
%            BBt(jj*nx+1:(jj+1)*nx , (ii-1)*nu+1:ii*nu) = AA*AiBt(:,:,ii);
%            AA = AA*A;
%         end   
%     end
%     mpc.BBt = BBt;
% 
%     C      = Cm;
%     Cc     = kron(eye(Nh),C);
%     
%     toc
% else
% %    code by Pablo
%     tic
%     L = zeros(nx*Nh,nx);
%     S = zeros(nx*Nh,nu*Nh);
%     Ai = Am;
%     Bi = Bm;
%     L(1:nx,1:nx) = Am;
%     S(1:nx,1:nu) = Bm;
%     for idx = 1:Nh-1
%         L(idx*nx+1:(idx+1)*nx,1:nx) = Ai*L((idx-1)*nx+1:idx*nx,1:nx);
%         S(idx*nx+1:(idx+1)*nx,1:(idx+1)*nu) = [Ai*S((idx-1)*nx+1:idx*nx,1:idx*nu),Bi];
%     end
%     C      = Cm;
%     Cc     = kron(eye(Nh),C);
%     toc
% end
% 
% mpc.AA  = L;
% mpc.BBt = S;
% mpc.CC  = Cc;

end

