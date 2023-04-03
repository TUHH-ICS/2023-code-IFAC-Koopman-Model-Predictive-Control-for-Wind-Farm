function mpc = MPCinit_ctrl(sol,Wp,mpc)
% MPCinit_ctrl initializes controller model params cf,cp, bcoef

mpc.cf           = sol.turbine.cf;
mpc.cp           = sol.turbine.cp;

for kk = 1:Wp.turbine.N    
    mpc.bcoef{kk} = mpc.num(2)*[-mpc.cf(kk);mpc.cp(kk);1]; %changing init with zeros                                   
end    

