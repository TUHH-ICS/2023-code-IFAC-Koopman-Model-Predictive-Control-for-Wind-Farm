function mpc = MPCinit_ctrl(sol,Wp,mpc)
% controller models                     
mpc.cf           = sol.turbine.cf;
mpc.cp           = sol.turbine.cp;

for kk = 1:Wp.turbine.N    
    mpc.bcoef{kk} = mpc.num(2)*[-mpc.cf(kk);mpc.cp(kk);1]; %changing init with zeros                                   
end    

