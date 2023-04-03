function [u,Koop,XX,Uf,lcp,xkoop,add,timing]  = KqLMPC2(Koop,x,x_,MPC,Uf,XX,lcp2,options)
%% Evaluate if dictionary has to be updated

timing = zeros(4,1);
start_timing = uint64(4);
start_timing(1) = tic();   % total time in function
% Tuning knobs
v = Koop.v;      % Maybe later v=0.2 if MPCenable????

% number of observables because of states and inputs
nx = options.nx;
ni = options.ni;
n = nx + ni;
nx_v = nx * 2;  % velocity based approach, twice number of states

% Get Koopman variables
A = Koop.A;
G = Koop.G;
K = Koop.K;
P = Koop.P;

% Get MPC varaiables
if isfield(MPC,'Q_')
Q_      = MPC.Q_;
R_      = MPC.R_;
C_delta = MPC.C_delta;
Ulim    = MPC.Ulim;
% Caug    = MPC.Caug;
% ylim    = MPC.ylim;
ref     = MPC.ref;
Uold    = kron(ones(options.N_Prediction_Horizon,1),MPC.Uold);
end

start_timing(2) = tic();  % time for updating library

% generate data for comparison
Measured_X_data = zeros(nx, options.Library_Update_Prediction_Horizon);    % initilize matrix for measured data
Predicted_X_from_data_ext = zeros(nx, options.Library_Update_Prediction_Horizon + 1);   % initilize matrix for predicted data

Predicted_X_from_data_ext(:,1) = x_(1:nx,1);    % Store as start point for prediction the last value before the horizon starts
Khat = K(:,1:nx)';  % Truncate the Koopman operator that it can be used for directly calculating the system output

% Fill the data matrices with measured and predicted data
for i = 0:options.Library_Update_Prediction_Horizon - 1
    Measured_X_data(:,i+1) = x(i*n+1:(i+1)*n-ni);   % Store further data from old to young
    Predicted_X_from_data_ext(:,i+2) = Khat * K_psi([Predicted_X_from_data_ext(:,i+1); x_((i+1)*nx+1 : (i+1)*nx+ni)]);    % Khat * Psi([xk_-1;uk_-1])
end

Predicted_X_from_data = Predicted_X_from_data_ext(:,2:end); %Remove first column from Matrix since it was only used for the initial prediction

% Check if inf-norm is larger the v of the signal horizon for any state
% sequence
error_norm_channels = Koop.error_norm_channels;
error_norm = norm(Measured_X_data(error_norm_channels,:) - Predicted_X_from_data(error_norm_channels,:),'inf');
% if norm(PHX([1 2 3 5 6 7])-KPHX([1 2 3 5 6 7]),'inf') > v
% error_norm = norm(Measured_X_data(error_norm_channels)-Predicted_X_from_data(error_norm_channels),'inf');
add = [error_norm;0];    % Output if data was added to the library
if error_norm > v
    P = P + 1;  % Number of library updates
    psi_xk_1 = K_psi(x_(end-n+1:end));   % predict xk-1
    psi_xk = K_psi(x(end-n+1:end));      % predict xk
    G = 1/P*((P-1)*G + psi_xk_1*psi_xk_1'); % Update G
    A = 1/P*((P-1)*A + psi_xk_1*psi_xk');   % Update A
    K = pinv(G)*A;  % Update Koopman operator
    add = [error_norm;1];   % Store the error_norm and that the Koopman operator was updated
    
    Koop.A = A;
    Koop.G = G;
    Koop.K = K;
    Koop.P = P;
end
Khat = K(:,1:nx)';  % Update the truncated Koopman operator

%  
ny = 2;
dt = 1;
noStates = 6;
approxA = K(1:noStates,1:noStates); % system matrix
approxB = K(1:noStates,noStates+1:end);
approxC = [eye(ny),zeros(ny,noStates-ny)];
approxD = zeros(ny,ni);
sysC = d2c(ss(approxA,approxB,approxC,approxD,dt),'tustin');
add(3) = max(real(eig(sysC)));

timing(2) = toc(start_timing(2));  % time for updating library

%% Do a new prediction with one step ago predicted u
NN = options.N_Prediction_Horizon;

% Initlize Matrices and first prediction
ph = K_psi(x(end-n+1:end)); % Calculate observables        
XXk = zeros(nx*NN,1);       % Storage Matrix for all predicted xk+i
XXk(1:nx,1) = Khat*ph;         % Calculate one step ahead prediction 

% further predictions
for i = 1:options.N_Prediction_Horizon - 1
    ph = K_psi([XXk((i-1)*nx+1:i*nx); Uf((i-1)*ni+1:i*ni)]);
    XXk(i*nx+1:(i+1)*nx) = Khat*ph;
end
xkoop = XXk(end-nx+1:end);  % output last prediction using Koopman

%% qLMPC
start_timing(3) = tic();  % Complete qLMPC algorithm
if options.MPC_enable
    %% MPC tuning knobs
    iter = options.qLMPC_iterations;
%     % Q_      = blkdiag(kron(eye(NN-1),diag([1 50 100 .1 1 1 1 0 0 0 0 0 0 0])),10*diag([1 50 100 .1 1 1 1 0 0 0 0 0 0 0])); %ulim 0.1
%     Q_      = blkdiag(kron(eye(NN-1),diag([1 100 100 .01 5 5 5 0 0 0 0 0 0 0])),10*diag([1 100 100 .01 5 5 5 0 0 0 0 0 0 0])); 
% %     R_      = kron(eye(NN),diag([5000,750]));
%     R_      = kron(eye(NN),10000*eye(ni));
%     ref     = kron(ones(NN,1),[0;r(1);r(2);45;0;0;0;0;0;0;0;0;0;0]);
%     C_delta = tril(kron(ones(NN),eye(ni)));
%     Uold    = kron(ones(NN,1),x(end-ni+1:end));
%     Ulim    = kron(ones(NN,1),[0.5;2]); % not maximum [0.66; 2.44]
%     Caug    = kron(eye(NN),[1 0 0 0 0 0 0 0 0 0 0 0 0 0]);
%     ylim    = pi/2;
    % Q_      = blkdiag(kron(eye(NN-1),diag([1 50 100 .1 1 1 1 0 0 0 0 0 0 0])),10*diag([1 50 100 .1 1 1 1 0 0 0 0 0 0 0])); %ulim 0.1

    
    %% Get previous data / Initilization of variables   
    u = x(end-ni+1:end);    % get previous u
    dU = zeros(ni*NN,1);       % initialize dU
    xk = x(end-n+1:end-ni); % get xk
    xk_ = x_(end-n+1:end-ni);   % get xk-1
    
    S = zeros(nx_v*NN,ni*NN);   % Initialize S matrix 
    H = zeros(nx_v*NN,nx_v);    % Initialize H matrix

    % qLMPC Scheme
    for it = 1:iter
        %% Preparation of Matrices
        Jx = dphi_dx([xk;u]);   % Calculation of derivatives
        Ju = dphi_du([xk;u]);   % Calcluation of derivatives

        Ak = [eye(nx)  Khat*Jx;  % Create A-Matrix of velocity approach
              zeros(nx) Khat*Jx];
        B = [Khat*Ju;Khat*Ju];  % Create B-Matrix of velocity approach
 
        H(1:nx_v,1:nx_v) = Ak;  % 1st step to fill the H matrix
        S(1:nx_v,1:ni) = B;     % 1st step to fill the S matrix

        % further step for filling the H and S matrices
        for mm=1:NN-1 
            Jx = dphi_dx([XX((mm-1)*nx_v+1:(mm-1)*nx_v+nx); Uf((mm-1)*ni+1:mm*ni)]); % Calculation of derviatives
            Ju = dphi_du([XX((mm-1)*nx_v+1:(mm-1)*nx_v+nx); Uf((mm-1)*ni+1:mm*ni)]); % Calculation of derviatives

            Ak = [eye(nx)  Khat*Jx;  % Create A-Matrix of velocity approach
                  zeros(nx) Khat*Jx];
            B = [Khat*Ju;Khat*Ju];  % Create B-Matrix of velocity approach

            H(mm*nx_v+1:(mm+1)*nx_v,1:nx_v) = Ak*H((mm-1)*nx_v+1:mm*nx_v,:);  % fill the H matrix
            S(mm*nx_v+1:(mm+1)*nx_v,1:mm*ni) = Ak*S((mm-1)*nx_v+1:mm*nx_v,1:mm*ni); 
            S(mm*nx_v+1:(mm+1)*nx_v,mm*ni+1:(mm+1)*ni) = B;             % fill the S matrix
        end

        sQ = S'*Q_;
        H_qp = 2*(sQ*S+R_);
        g = 2*sQ*(H*[xk;xk-xk_]-ref);
        
        % Calculate Constraint matrices
        % Ac = [-C_delta;C_delta;-[Caug;-Caug]*S]; ??????
        % bc = [Ulim+Uold;Ulim-Uold;-[-ylim-Caug*L*[xk;xk-xk_];-ylim+Caug*L*[xk;xk-xk_]]];

        Ac = [-C_delta;C_delta];
        bc = [Ulim+Uold;Ulim-Uold];
        
        timing(3) = toc(start_timing(3));  % Create matrices
        
        lcp = lcp2;
        start_timing(4) = tic();  % Just solver time
        % Selection of solver
        switch options.solver
            case 0 % quadprog
                dU = quadprog(H_qp,g,Ac,bc);      
            case 1 % qpOASES
                dU = qpOASES(H_qp,g,Ac,[],[],[],bc);%,opt);
            case 2 % lemke
                Ac_ = -Ac;
                bc_ = -bc;
                Hi = eye(size(H_qp,1))/H_qp;
                M = Ac_*Hi*Ac_';
                q = -Ac_*Hi*g-bc_;
                lcp = lemke(M,q,lcp2);
                dU = Hi*(-Ac'*lcp-g);
            case 3 % QPsol
                dU = QPsol(H_qp,g,[],[],Ac,bc);
        end
        timing(4) = toc(start_timing(4)); % Just solver time
        
        dU(isnan(dU)) = 0;  % Set NAN to 0
        XX = H*[xk;xk-xk_] + S*dU;  % Calculate predictions of velocity states using qLMPC and most recent calculated U

        u = u + dU(1:2); % Update u for the next step
        Uf = C_delta * dU + Uold;
    end
else
    % Calculate nesseary variables in the learning phase
    u = zeros(ni,1);    % Set controller output to 0
    lcp = lcp2;         % Forward the previous lcp solution
    XX(1:nx_v) = [XXk(1:nx);XXk(1:nx)-x(end-n+1:end-ni)];   % Calculate predictions of velocity states based on Koopman predictions and actuall measurements
    for i = 1:NN-1
        XX(i*nx_v+1:(i+1)*nx_v) = [XXk(i*nx+1:i*nx+nx);XXk(i*nx+1:i*nx+nx)-XXk((i-1)*nx+1:(i-1)*nx+nx)];
    end
end

timing(1) = toc(start_timing(1)); % total time in function