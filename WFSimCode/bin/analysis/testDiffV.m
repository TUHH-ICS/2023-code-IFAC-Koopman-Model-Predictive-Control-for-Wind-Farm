 KWfsim = Wp.turbine.powerscale * 2* 0.5*rho* A; % Wp.turbine.powerscale = 0.9500
 kGauss = 1/2* rho * A* Cp(1);
 
 Ur2Plot =  [3.4199    3.5347    3.7231    3.9295    4.1300    4.3033    4.4538    4.5888]

 
 K2K = Wp.turbine.powerscale * 2/Cp(2);
 
 Utilde = (K2K * Ur2Plot.^3) .^ (1/3)

 Utilde = (K2K * Ur2Plot.^3) .^ (1/3)

Utilde =

    5.0460    5.2155    5.4935    5.7980    6.0939    6.3495    6.5716    6.7708
     0.6308    0.6519    0.6867    0.7248    0.7617    0.7937    0.8215    0.8464
     
     
     
   