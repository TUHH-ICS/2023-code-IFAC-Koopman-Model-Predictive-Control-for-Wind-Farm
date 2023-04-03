function turbInputSet = controlSet_ref12(Wp)
    turbInputSet = struct();
    
    nTurbs = length(Wp.turbine.Crx);
    turbInputSet.t = 0:Wp.sim.h:20e2; % 10,000 seconds of simulation
    NN = (length(turbInputSet.t)-1)/2;
    turbInputSet.phi = zeros(nTurbs,length(turbInputSet.t));
    turbInputSet.CT_prime = [0.2*ones(nTurbs-1,NN/2.5),.5*ones(nTurbs-1,NN/2.5),1*ones(nTurbs-1,NN/2.5),1.5*ones(nTurbs-1,NN/2.5),2*ones(nTurbs-1,NN/2.5+1);...
        0.2*ones(nTurbs-1,length(turbInputSet.t))];
    %turbInputSet.CT_prime(2,length(turbInputSet.t)) = 2.0*ones(nTurbs-1,length(turbInputSet.t));
    turbInputSet.interpMethod = 'nearest'; % Nearest since no interpolation necessary
end

