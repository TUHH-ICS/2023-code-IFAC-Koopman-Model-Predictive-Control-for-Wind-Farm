function [Wp,sol,sys] = InitWFSim(Wp,options,plotMesh)

%% INITWFSIM  Initializes the WFSim model

    %Input%%
    % Wp : Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list 
    % options: Choose model solver options
    % plotMesh : Plot mesh, turbine locations, and print grid offset values. Default: false.
    %Output%%
    % Wp : Smulation layout afer updating the turbine loctaion and meshing 
    % sol : Contains information of Flow field(unform), flow velocity and pressure
    % sys : 
    
    
%% 
    % Create empty structs
    sys = struct; % This struct will contain all the system matrices at time k
       
    % Import simulation scenario (meshing, atmospheric properties, turbine settings)
    [Wp] = meshing(Wp,plotMesh,plotMesh); 

    % Initialize time vector for sol at time k = 0
    sol = struct('k',0,'time',0);% This struct will contain the solution (flowfields, power, ...) at time k
    
    % Initialize flow fields as uniform ('no turbines present yet')nw each
    % grid as the value of velocity(x&y direction) and pressure.
    [sol.u,sol.uu] = deal( Wp.site.u_Inf(1) *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );% Y = deal(X) copies the single input to all the requested outputs.  
    [sol.v,sol.vv] = deal( Wp.site.v_Inf(1) *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.p,sol.pp] = deal( Wp.site.p_init * ones(Wp.mesh.Nx,Wp.mesh.Ny) ); 

    % Initialize the linearized solution variables, if necessary
    if options.Linearversion
        sol.ul = sol.u;
        sol.vl = sol.v;
        sol.pl = sol.p;
        [sol.du,sol.dv,sol.dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
    end

    % Compute boundary conditions and system matrices B1, B2.
    [sys.B1,sys.B2,sys.bc]  = Compute_B1_B2_bc(Wp);
    sys.pRCM                = []; % Load empty vector
    
    % Compute projection matrices Qsp and Bsp. These are only necessary if
    % the continuity equation is projected away (2015 ACC paper, Boersma).
    % todo wat he need to do
    if options.Projection
        [sys.Qsp, sys.Bsp]  = Solution_space(sys.B1,sys.B2,sys.bc);
        Wp.Nalpha           = size(sys.Qsp,2);
    end
end

