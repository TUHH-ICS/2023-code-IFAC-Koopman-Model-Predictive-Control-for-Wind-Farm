function [output,sol] = Actuator(Wp,sol,options)
%% Import variables
Nx              = Wp.mesh.Nx;
Ny              = Wp.mesh.Ny;
dyy2            = Wp.mesh.dyy2;
xline           = Wp.mesh.xline;
yline           = Wp.mesh.yline;
ylinev          = Wp.mesh.ylinev;
Rho             = Wp.site.Rho;
Drotor          = Wp.turbine.Drotor;
powerscale      = Wp.turbine.powerscale;
N               = Wp.turbine.N;
F               = Wp.turbine.forcescale; 
input           = sol.turbInput;
Projection      = options.Projection;
Linearversion   = options.Linearversion;

%% Scale power and force with yaw angel
% Comparing Three Aerodynamic Models for Predicting The Thrust and Power Characteristics of Yawed Floating Wind Turbine Rotors
phiAbs = abs(sol.turbInput.phi);
ctVec = [0,15,30,45;0.79, 0.76,0.69,0.58]; %temp = interp1([0,15,30,45],[0.79, 0.76,0.69,0.58],0:45);
ctVec(2,:) = ctVec(2,:)/ctVec(2,1);
F =  F.* interp1(ctVec(1,:),ctVec(2,:),phiAbs(:));

cpVec = [0,15,30,45; 0.5, 0.48,0.4,0.26]; %temp = interp1([0,15,30,45],[0.5, 0.48,0.4,0.26],0:45);
cpVec (2,:) = cpVec(2,:)/cpVec(2,1);
powerscale =  powerscale.* interp1(cpVec(1,:),cpVec(2,:),phiAbs(:));

% F  = repmat(Wp.turbine.forcescale,N,1);
% powerscale = repmat(Wp.turbine.powerscale,N,1);

%%
Ar              = pi*(0.5*Drotor)^2;

[Sm.x,Sm.dx]    = deal(sparse(Nx-3,Ny-2));            % Input x-mom nonlinear and linear
[Sm.y,Sm.dy]    = deal(sparse(Nx-2,Ny-3));            % Input y-mom nonlinear and linear
[Sm.xx,Sm.dxx]  = deal(sparse((Nx-3)*(Ny-2),2*N));    % Input x-mom nonlinear and linear qlpv
[Sm.yy,Sm.dyy]  = deal(sparse((Nx-2)*(Ny-3),2*N));    % Input y-mom nonlinear and linear qlpv


if Linearversion
    Smdu            = sparse(Nx-3,Ny-2);
    Smdv            = sparse(Nx-2,Ny-3);
end

for kk=1:N
   
    x  = xline(kk,:);  % Turbine x-pos in field
    y  = yline{kk};    % Turbine y-pos in field
    % yv = ylinev{kk};   % Corrected turbine y-pos in field
    
    vv            = 0.5*(sol.v(x,y+1)-sol.v(x,y))+sol.v(x,y); %0.5*diff(sol.v(x,yv))+sol.v(x,yv(1:end-1)); 
    uu            = sol.u(x,y);
    U{kk}         = sqrt(uu.^2+vv.^2);
    U_v(kk)       = mean(U{kk});
    phi{kk}       = atan(sol.v(1,1)/sol.u(1,1));
    
    phiCell =     atan(vv./uu); 
    Ue{kk}        = cos(input.phi(kk)/180*pi + phiCell).*U{kk};
    Ur(kk)        = mean(Ue{kk});  
    CT(kk)        = input.CT_prime(kk);                     % Import CT_prime from inputData
    Phi(kk)       = input.phi(kk);
    
    %% Thrust force     
    Fthrust         = F(kk)*1/2*Rho*Ue{kk}.^2*CT(kk);           % Using CT_prime
    Fx              = Fthrust.*cos(phi{kk}+Phi(kk)*pi/180);
    Fy              = Fthrust.*sin(phi{kk}+Phi(kk)*pi/180);
    cf(kk)          = F(kk)*1/2*Rho.*mean(dyy2(1,y));
    Force(kk)       = -mean(Fx.*dyy2(1,y));
    
    %% Power
    cp(kk)          = powerscale(kk)*.5*Rho*Ar;
    Power(kk)       = cp(kk)*CT(kk)*mean(Ue{kk}.^3);    
    
    %% Input to Ax=b
    Sm.x(x-2,y-1)           = -Fx'.*dyy2(1,y)';               % Input x-mom nonlinear                           
    Sm.y(x-1,y(2:end)-2)    = Fy(2:end)'.*dyy2(1,y(2:end))';  % Input y-mom nonlinear
    
    % Apply the force to the trailing cells to achieve a higher (LES-like) wake deflection
    Sm.y(x,y(2:end)-2)   = Fy(2:end)'.*dyy2(1,y(2:end))';  
    Sm.y(x+1,y(2:end)-2) = Fy(2:end)'.*dyy2(1,y(2:end))';
    
    % Matrices for linear version
    if Linearversion
        dCT(kk) = input.dCT_prime(kk); 
        
        dFthrustdCT             = F*1/2*Rho*Ue{kk}.^2;
        dFxdCT                  = dFthrustdCT.*cos(phi{kk}+Phi(kk)*pi/180);
        dFydCT                  = dFthrustdCT.*sin(phi{kk}+Phi(kk)*pi/180);
         
        Sm.dx(x-2,y-1)          = -dFxdCT'*dCT(kk).*dyy2(1,y)';
        Sm.dy(x-1,y(2:end)-2)   = dFydCT(2:end)'*dCT(kk).*dyy2(1,y(2:end))';  
               
        dFdu                    = F*Rho*cos(Phi(kk)*pi/180)^2*CT(kk)*uu;
        dFdv                    = F*Rho*cos(Phi(kk)*pi/180)^2*CT(kk)*vv;     
        Smdu(x-2,y-1)           = -dFdu'.*dyy2(1,y)';
        Smdv(x-1,y-2)           =  dFdv'.*dyy2(1,y)';
        
        dSm.dx                  = blkdiag(diag(vec(Smdu')'),diag(vec(Smdv')'),sparse((Ny-2)*(Nx-2),(Ny-2)*(Nx-2)));
        
        % following for projection
        tempdx                  = sparse(Nx-3,Ny-2);
        tempdy                  = sparse(Nx-2,Ny-3);
        tempdx(x-2,y-1)         = -dFxdCT'.*dyy2(1,y)';
        Sm.dxx(:,kk)            =  vec(tempdx');                  % Input matrix (beta) x-mom linear       
        tempdy(x-1,y(2:end)-2)  = dFydCT(2:end)'.*dyy2(1,y(2:end))';
        Sm.dyy(:,kk)            = vec(tempdy');                   % Input (beta) y-mom linear qlpv
    end
    
    if Projection
        %% Input to qLPV
        % Clear for multiple turbine case
        tempx                   = sparse(Nx-3,Ny-2);
        tempy                   = sparse(Nx-2,Ny-3);
        tempx(x-2,y-1)          = -Fx'.*dyy2(1,y)';
        Sm.xx(:,kk)             = vec(tempx')/CT(kk);
        Sm.xx(:,N+kk)           = vec(tempx');
        if input.phi(kk)~=0
            Sm.xx(:,kk)         = Sm.xx(:,kk)/2;
            Sm.xx(:,N+kk)       = Sm.xx(:,N+kk)/(2*Phi(kk));         % Input x-mom nonlinear qlpv
        end
        
        tempy(x-1,y(2:end)-2)   = Fy(2:end)'.*dyy2(1,y(2:end))';
        Sm.yy(:,kk)             = vec(tempy')/CT(kk);
        Sm.yy(:,N+kk)           = vec(tempy');
        if input.phi(kk)~=0
            Sm.yy(:,kk)         = Sm.yy(:,kk)/2;
            Sm.yy(:,N+kk)       = Sm.yy(:,N+kk)/(2*Phi(kk));         % Input y-mom nonlinear qlpv
        end
        
     end
     
end

%% Write to outputs
sol.turbine.cp(:,1)       = cp;
sol.turbine.power(:,1)    = Power;
sol.turbine.CT_prime(:,1) = CT;
sol.turbine.Phi(:,1)      = Phi;
sol.turbine.Ur(:,1)       = Ur;
sol.turbine.U(:,1)       = U_v;
sol.turbine.cf(:,1)       = cf;
sol.turbine.force(:,1)    = Force;

output.Sm  = Sm;
if Linearversion>0
    output.dSm = dSm;
end