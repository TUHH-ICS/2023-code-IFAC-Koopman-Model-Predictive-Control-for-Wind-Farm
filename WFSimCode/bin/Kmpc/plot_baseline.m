function [mpc] = plot_baseline(Wp,sol,mpc)
%PLOT_BASELINE Summary of this function goes here
%   Detailed explanation goes here


delta_U = sol.turbInput.dCT_prime;
mpc.vect_R_baseline(sol.k) = delta_U' * mpc.R(1:9,1:9) * delta_U;

for nl = 1:Wp.turbine.N
    for nnl = 1:4
        velocity = sol.u(Wp.mesh.xline(nl),Wp.mesh.yline(nl){1});
    end
end

P_meas = ((pi*(Wp.turbine.Drotor)^2)/8)*(velocity*cos(turbInput.phi))^3*turbInput.CT_prime;
end

