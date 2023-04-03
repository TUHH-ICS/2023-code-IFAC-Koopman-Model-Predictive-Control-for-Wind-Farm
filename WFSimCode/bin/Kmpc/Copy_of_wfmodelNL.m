function [K, mpc] = wfmodelNL(sol,Wp,mpc,ll,K)


% for a wind turbine model state space representation.
% x(k+1) = Ax+b(Ct_prime)
% y=cx
% x = [Fk Pk C_hat_T_prime]'
% B = [C_t 0 0; 0 C_p 0; 0 0 1]
% C = [1 0 0: 0 1 0; 0 0 1]
% build B matrix of the turbine models for current times step

ny = 2; % real states Ur1 and Ur2
% if measured wind speed is used:

nTx = round(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx)-1;
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);
u1 = sol.u(nTx,nTy);

phiCos = cos(sol.turbine.Phi/180 *pi);

if size(K.B,2) == 3
    uIn = [sol.turbine.CT_prime; u1];
else
    uIn = sol.turbine.CT_prime;
end

nKoop = size(K.A,1);
polyKoop = 0; %K.PolyLiftingFunction;

mpc.structPrev.Ur1_prev1 = sol.turbine.Ur(1);
mpc.structPrev.Ur2_prev1 = sol.turbine.Ur(2);

if sol.k == 1
    
    mpc.structPrev.dUr1_prev1 = 0;
    mpc.structPrev.dUr2_prev1 = 0;
    mpc.structPrev.M1 = mpc.structPrev.Ur1_prev1; % Moving mean
    mpc.structPrev.M2 =mpc.structPrev.Ur2_prev1;
    mpc.structPrev.k = 1;
    [temp,~, mpc.structPrev]= K_psi([sol.turbine.Ur; uIn],polyKoop,nKoop,mpc.structPrev);
    K.xprev = temp(1:nKoop);
    mpc.xprev = K.xprev;
    
end

%K.Xf = zeros(length(mpc.xprev),mpc.Nh);
K.Xf1 = zeros(length(mpc.xprev),mpc.Nh);

for nn = 1:mpc.Nh
    if nn == 1
        K.Xf(:,nn) =  K.A*mpc.xprev + K.B*uIn;
        structPrev = mpc.structPrev;
       
        K.Xf1(:,nn) =  K.A*mpc.xprev+K.B*uIn;
        [tmp,~,structPrev] = K_psi([K.Xf1(1:ny,1); uIn],polyKoop,nKoop,structPrev);
        K.Xf1(:,nn) = tmp(1:size(K.A,1));
        mpc.structPrev = structPrev;
    else
        %K.Xf(:,nn) =  K.A*K.Xf(:,nn-1)+K.B*uIn;
        K.Xf1(:,nn) =  K.A*K.Xf1(:,nn-1)+K.B*uIn;
        [tmp,~,structPrev] = K_psi([K.Xf1(1:ny,nn); uIn],polyKoop,nKoop,structPrev);
        K.Xf1(:,nn) = tmp(1:size(K.A,1));
        
    end
end
%K.Xf =  K.A*mpc.xprev+K.B*sol.turbine.CT_prime;
%K.V = mpc.xprev(1:2,:);
K.V        = K.Xf1(1:2,:);for kk = 1:Wp.turbine.N
    
    if ll==1
        %K.V         = K.Xf(1:2,1);
        v           = repmat(sol.turbine.Ur(kk),mpc.Nh,1);%repmat(sol.turbine.Ur(kk),mpc.Nh,1);
        vest        = repmat(K.V(kk,1),mpc.Nh,1);
    else
        
        v          = repmat(sol.turbine.Ur(kk),mpc.Nh,1); %mpc.V(kk,:);%
        vest       = K.V(kk,:);
    end
    mpc.b{kk}   = [mpc.bcoef{kk}(1)*(v(1)* phiCos(kk))^2 mpc.bcoef{kk}(2)*(v(1)* phiCos(kk))^3 mpc.bcoef{kk}(3)]';
    mpc.best{kk}   = [mpc.bcoef{kk}(1)*(vest(1)* phiCos(kk))^2 mpc.bcoef{kk}(2)*(vest(1)* phiCos(kk))^3 mpc.bcoef{kk}(3)]';
    for nn = 1:mpc.Nh
        mpc.bt(:,kk,nn) = [mpc.bcoef{kk}(1)*(v(nn)* phiCos(kk))^2 mpc.bcoef{kk}(2)*(v(nn)* phiCos(kk))^3 mpc.bcoef{kk}(3)]';
        mpc.btest(:,kk,nn) = [mpc.bcoef{kk}(1)*(vest(nn)* phiCos(kk))^2 mpc.bcoef{kk}(2)*(vest(nn)* phiCos(kk))^3 mpc.bcoef{kk}(3)]';
    end
    
end

% K.Xf =  K.A*mpc.xprev+K.B*sol.turbine.CT_prime;
% K.xprev = K.Xf(:,1);
% build wind farm model
mpc.A = blkdiag(mpc.a{:});
mpc.B = blkdiag(mpc.b{:});
mpc.Best = blkdiag(mpc.best{:});
for kk = 1:mpc.Nh
    bCell  = cell(Wp.turbine.N,1);
    for nn = 1:Wp.turbine.N
        b{nn}     = mpc.bt(:,nn,kk);
        best{nn}  = mpc.btest(:,nn,kk);
    end
    mpc.Bt(:,:,kk) = blkdiag(b{:});
    mpc.Btest(:,:,kk) = blkdiag(best{:});
end
mpc.C = blkdiag(mpc.c{:});
mpc.xprev = K.Xf1(:,1);
end

