function [B1,B2,bc] = Compute_B1_B2_bc(Wp)

Nx     = Wp.mesh.Nx;% no of cell in x direction
Ny     = Wp.mesh.Ny;% no of cell in y direction
dxx2   = Wp.mesh.dxx2;%calculate the difference bw two secondary point in x direction(i) (ie e.g. ?x4,5) and sized to 100* 42 matrix.
dyy2   = Wp.mesh.dyy2;%%calculate the difference bw two secondary point in x direction(j) (ie e.g. ?x4,5) and sized to 100* 42 matrix.
Rho    = 1; % Density of  air
%Rho    = 1.225;%Wp.site.Rho;

u_Inf  = Wp.site.u_Inf(1); % Init. Longitudinal Wind speed. (m/s)

% system matrix in B1 and B2 in eq 20 in matrix Ek matrix
Bm1                  = Rho*(spdiags(-ones((Nx-2)*(Ny-2),1).*vec(dyy2(2:end-1,2:end-1)'),0,(Nx-3)*(Ny-2),(Nx-2)*(Ny-2))+...
    spdiags(ones((Nx-2)*(Ny-2),1).*vec(dyy2(2:end-1,2:end-1)'),Ny-2,(Nx-3)*(Ny-2),(Nx-2)*(Ny-2)));
Bm2                  = Rho*(spdiags(-ones((Nx-2)*(Ny-2),1).*vec(dxx2(2:end-1,2:end-1)'),0,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2))+...
    spdiags(ones((Nx-2)*(Ny-2),1).*vec(dxx2(2:end-1,2:end-1)'),1,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2)));
Bm2(Ny-2:Ny-2:end,:) = [];

B1 = Bm1'; B2 = Bm2';

bc          = zeros((Ny-2)*(Nx-2),1);
bc(1:Ny-2)  = -Rho*u_Inf*dyy2(1,2:end-1)';% todo unit of dyy is m

B1((Ny-2)*(Nx-3)+1:(Ny-2)*(Nx-2),:) = 0; % u_{Nx,J}=u_{Nx-1,J}

for kk=1:Nx-2
    B2(kk*(Ny-2),:)= 0; % v_{I,Ny}=v_{I,Ny-1}
end

for kk=0:Nx-3
    B2(kk*(Ny-2)+1,:)= 0; % v_{I,3}=v_{I,2} for I=2,3,...,Nx
end

B1 = B1';
B2 = B2';




