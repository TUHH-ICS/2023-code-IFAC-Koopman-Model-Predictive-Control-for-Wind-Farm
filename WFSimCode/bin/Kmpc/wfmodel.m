function [K, mpc] = wfmodel(sol,Wp,mpc,ll,K)

% for a wind turbine model state space representation.
% x(k+1) = Ax+b(Ct_prime)
% y=cx
% x = [Fk Pk C_hat_T_prime]'
% B = [C_t 0 0; 0 C_p 0; 0 0 1]
% C = [1 0 0: 0 1 0; 0 0 1]
% build B matrix of the turbine models for current times step

K.Xf = zeros(length(mpc.xprev),mpc.Nh);
for nn = 1:mpc.Nh
    if nn == 1
        K.Xf(:,nn) =  K.A*mpc.xprev+K.B*sol.turbine.CT_prime;
    else
        K.Xf(:,nn) =  K.A*K.Xf(:,nn-1)+K.B*sol.turbine.CT_prime;
    end
end
%K.Xf =  K.A*mpc.xprev+K.B*sol.turbine.CT_prime;
%K.V = mpc.xprev(1:2,:);
K.V        = K.Xf(1:2,:);
for kk = 1:Wp.turbine.N
    
    if ll==1
        %K.V         = K.Xf(1:2,1);
        v           = repmat(sol.turbine.Ur(kk),mpc.Nh,1);%repmat(sol.turbine.Ur(kk),mpc.Nh,1);
        vest        = repmat(K.V(kk,1),mpc.Nh,1);
    else
        
        v          = repmat(sol.turbine.Ur(kk),mpc.Nh,1); %mpc.V(kk,:);%
        vest       = K.V(kk,:);
    end
    mpc.b{kk}   = [mpc.bcoef{kk}(1)*v(1)^2 mpc.bcoef{kk}(2)*v(1)^3 mpc.bcoef{kk}(3)]';
    mpc.best{kk}   = [mpc.bcoef{kk}(1)*vest(1)^2 mpc.bcoef{kk}(2)*vest(1)^3 mpc.bcoef{kk}(3)]';
    for nn = 1:mpc.Nh
        mpc.bt(:,kk,nn) = [mpc.bcoef{kk}(1)*v(nn)^2 mpc.bcoef{kk}(2)*v(nn)^3 mpc.bcoef{kk}(3)]';
        mpc.btest(:,kk,nn) = [mpc.bcoef{kk}(1)*vest(nn)^2 mpc.bcoef{kk}(2)*vest(nn)^3 mpc.bcoef{kk}(3)]';
    end
    
end

%K.Xf =  K.A*mpc.xprev+K.B*sol.turbine.CT_prime;
K.xprev = K.Xf(:,1);
% build wind farm model
mpc.A = blkdiag(mpc.a{:});
mpc.B = blkdiag(mpc.b{:});
mpc.Best = blkdiag(mpc.best{:});
for kk = 1:mpc.Nh
    for nn = 1:Wp.turbine.N
        b{nn}     = mpc.bt(:,nn,kk);
        best{nn}     = mpc.btest(:,nn,kk);
    end
    mpc.Bt(:,:,kk) = blkdiag(b{:});
    mpc.Btest(:,:,kk) = blkdiag(best{:});
end
mpc.C = blkdiag(mpc.c{:});
mpc.xprev = K.xprev;
end

