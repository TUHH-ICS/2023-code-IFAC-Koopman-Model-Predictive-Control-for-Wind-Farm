function [mpc] = plot_qLPVMPC(mpc,sol,nlIdx)
%PLOT_QLPVMPC Summary of this function goes here
%   Detailed explanation goes here

mpc.vect_Q_qlpvmpc(sol.k,nlIdx) = mpc.vect_E'*mpc.Q(1,1)*mpc.vect_E;
mpc.vect_R_qlpvmpc(sol.k,nlIdx) =  mpc.vect_dU'*mpc.R(1:9,1:9)* mpc.vect_dU;



end

