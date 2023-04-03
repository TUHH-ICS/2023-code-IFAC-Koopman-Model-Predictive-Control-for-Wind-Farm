function [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd,x] = dynamicmodedecomposition(states, Inputs, Outputs, Deterministic,method,r,maindir)

% This function aims to build a reduced order model from the states,
% input/output information and deterministic states gathered in the
% simulation and resampled

%% INPUT ARGUMENTS
%states: states to build DMD matrices
%inputs: [U1 U2 U3 ... Un];
%outputs: [Y1 Y3 .... Yn];
%deterministic: [Xd1 Xd2 Xd3 Xd4 ... Xdn];
%method: which method to use for DMD
%r: truncation level (number of singular values to retain)


%% OUTPUTS ARGUMENTS
% sys_red: the state space systems, according to the number of singular
% values used
% FITje: the Fit of the model for the given data
% U, S, V: the matrices resulting from SVD to be used for DMD states reconstruction
%method: method is later used, as reconstruction depends on the method
%used, even though methodology is similar
% X, X_p: DMD matrices to be used later for reconstruction

%% DMD - Dynamic Mode Decomposition

% x(k+1) ~ A x(k)
% X' ~ AX

% X = [   |  |        |  ]
%     [   x1 x2 ... xm-1 ]
%     [   |  |        |  ]

% X'= [   |  |        |  ]
%     [   x2 x3 ...  xm  ]
%     [   |  |        |  ]

%define necessary matrices for DMD
X      = states(:,1:end-1);
X_p    = states(:,2:end);

out      = Outputs(:,1:end-1);
inp      =Inputs(:,1:end-1);

Xd     = Deterministic(:,1:end-1);
Xd_p   = Deterministic(:,2:end);

%% (0) DMD - there is no control action
% part = 2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);

approxA = cell(r,1);
approxB = cell(r,1);
approxC = cell(r,1);
approxD = cell(r,1);
sys_red = cell(r,1);
 

if method==0 %dmd algortihm
    
    dirdmd = 'DMDresults_DMD';
    dirdmd = fullfile(maindir,dirdmd);
    if ~exist(dirdmd,'dir')
        mkdir(dirdmd);
    end
    
    [U, S, V] = svds(X,r);
    
    for si=1:1:r
        
        % part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);
        
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
        
        approxA{si} = Util'*(X_p)*Vtil/Sigtil; % *inv(Sigtil);
        approxB{si} = zeros(si,1);
        approxC{si} = zeros(2,si);
        approxD{si} = zeros(2,1);
        sys_red{si} = ss(approxA{si},approxB{si},approxC{si},approxD{si},2);
        
        FITje=0;
    end
    
    
elseif method==1 %dmdc algortihm
    %%  (1) DMD - there is a extrnal forcing term
    
    dirdmd='DMDresults_DMDc';
    dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir')
        mkdir(dirdmd);
    end
    
    %Goal of DMDc is to characterize the relationship between three
    %measurments: current measurment x(k), the future state x(k+1) and the
    %current control u(k). Relationship is approximated by the canonical
    %discrete linear dynamical system.
    %Assumptions: all states are observable
    %are C will be identity.
    
    % x(k+1) ~  Ax(k) + Bu(k)
    % y(k)   ~  Cx(k) + Du(k)
    
    % X' ~ AX + BU
    % Y  ~ CX + DU
    
    % X' ~ AX + BU ~ [A B][ X ] = G?
    %                     [ U ]
    
    %2 Perform SVD on the augmented data matrix SVD(?)=USV such that
    %G=X'VS^(-1)U*
    
    %3. Separate A and B by splitting left singular vectors into 2 separate
    %components
    
    % [A B] ~ [X'VS^(-1)U*1, X'VS^(-1)U*2 ]
    
    % U*1:
    % U*2:
    
    %4. Construct reduced order subspace from the output measurments X'. U
    %cannot be used to find low rank model of the dynamics and input matrixes
    %since it is defined for the input space which now includes both the state
    %measurments and the exogeneous inputs
    
    %SVD(X')=?S^V*
    
    % ?  = ?*A? = ?*X'VS^(-1)U*1 ?
    % B~ = ?*B  = ?*X'VS^(-1)U*2
    
    Omega = [X;inp];
    
    [U, Sig, V] = svds(Omega,r);
    [Uf, ~, ~] = svds(X_p,r);
    %     FITje = zeros(2, r);
    %     OMEGA = {};
    %     DAMPING ={};
    %     x=cell(r,1);
    
    for si = 1:r
        
        % part = 2; subpart = 5; f = MPC_progress(part,subpart,f,si,r);
        
        Util=U(:,1:si);
        Sigtil=Sig(1:si,1:si);
        Vtil=V(:,1:si);
        
        Uhat=Uf(:,1:si);
        % Sighat=Sigf(1:si,1:si);
        % Vbar=Vf(:,1:si);
        
        n=size(X,1);
        q=size(inp,1);
        U_1=Util(1:n,:);
        U_2=Util(n+1:n+q,:);
        
        approxA{si} = Uhat'*(X_p)*Vtil/Sigtil*U_1'*Uhat;
        approxB{si} = Uhat'*(X_p)*Vtil/Sigtil*U_2';
        approxC{si} = eye(si,si);
        approxD{si} = zeros(si,1);
        sys_red{si} = ss(approxA{si},approxB{si},approxC{si},approxD{si},2);
        
        
        close all
    end
    xo = dinit(sys_red{r}.A,sys_red{r}.B,sys_red{r}.C,sys_red{r}.D,Inputs',(Uhat'*states)');
    ysim = lsim(sys_red{r}, Inputs',[],xo);
    x=ysim;
    Xd={};
    FITje=0;
    
    
    %% (2) ioDMD: Input Output Dynamic Mode Decomposition
    
elseif method == 2 %ioDMD
    
    %The goal of ioDMD is to capitalize on DMDc and extend it so a full
    %state space system may be obtained via usual subspace system
    %identification methods, estimating matrices A,B,C,D via lieast squares
    
    dirdmd = 'DMDresults_IODMD';
    dirdmd = fullfile(maindir,dirdmd);%change from strcut
    if ~isfolder(dirdmd)
        mkdir(dirdmd);
    end
    
    % TOdOtest SVD
%         [~,Sf] = svds(double(X),908);%[~,Sf]=svd(X');
%         diagSf = diag(Sf);
%         Sfint = cumsum(diagSf);
%         Sfsum = sum(diagSf);
%         SfintNorm=Sfint/Sfsum;
%         figure; plot(SfintNorm) % SfintNorm(100) is 100 value is 80%
%         vec95 = min(find(SfintNorm>0.95));% 396 should be taken
%         vec90 = min(find(SfintNorm>0.90));%244
%         vec99 = min(find(SfintNorm>0.99));% 701
    [U,S,V] = svds(X,r);% TODo plot all singulatr values
    FITje=zeros(2, r);
    %FITje_fil=zeros(2, r);
%     OMEGA={};
%     DAMPING={};
    x = cell(r,1);
%     [U,S,V] = svd(X);
     %X1 = U *S *V';
%     
%     Vt = V';
%     VforRelieff = Vt(1:size(X,1),:);
    
    %[idx,weights] = relieff(X',out(2,:)',10);
    
    %[idx,weights] = relieff( X',out(2,:)',10);

     %V=V(idx,:);
     %X_p=X_p(idx,:);
%     [U,S,V] = svds(X(idx(1:100),:),r);
%     X_p = X(idx(1:100),:);
%     states = X_p;

    for si = 1: r 
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
        % [A^~ B^~; C^~ D^~]=[U^{-1}X(k+1);y(k)][sigmaV^{T};Gamma(k)](7)
        all=[Util'*X_p;out]*pinv([Sigtil*Vtil';inp]);
        % A^~ B^~; C^~ D^~
        si1 = si;
        approxA{si1}=all(1:si,1:si);
        approxB{si1}=all(1:si,si+1:end);
        approxC{si1}=all(si+1:end, 1:si);
        approxD{si1}=all(si+1:end, si+1:end);       
        sys_red{si1}=ss(approxA{si1},approxB{si1},approxC{si1},approxD{si1},2);%TODO Ts should be in main
        
        [FITje,fig1,x] = evaluatemodel(sys_red,si1,Inputs,Outputs,FITje,'identification',x,states,U,Deterministic,method,0);
        
        if ~isempty(fig1)
        warning off
        export_fig(fig1,strcat(dirdmd,'/image',num2str(10000+si1)),'-nocrop','-m2')
        warning on
        close all
        end
    end
    
    fig200 = VAFpermodes(FITje,r,{});
    warning off
    export_fig(fig200,strcat(dirdmd,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
    warning on
    Xd={};
    
    %% (3) extioDMD: Extended Input Output Dynamic Mode Decomposition
    
elseif method == 3
    
    % The goal of extioDMD is to obtain a state space system, as ioDMD
    % performs, but etending the existing state to others which are not
    % necessarily related with the previous. This cpaitalizes on the
    % convergence of DMD results and the Koopman Operator, where there is
    % significant evidence that by incuding non linear observables (which
    % may be functions of the pre exisitng states, or not) DMD provides
    % better results.
    %These observables used to extend the current state space are referred
    %to as determinisitc states, as they are measurable and known for the
    %current scenario
    
    dirdmd='DMDresults_EIODMD';
    dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir')
        mkdir(dirdmd);
    end
    
    dirdmdident='DMDresults_EIODMD/ident';
    dirdmdident=strcat(maindir,dirdmdident);
    if ~exist(dirdmdident,'dir')
        mkdir(dirdmdident);
    end
    
    [U,S,V]=svds(X,r);
    FITje=zeros(2, r);
    % FITje_fil=zeros(2, r);
%     OMEGA={};
%     DAMPING={};
    x=cell(r,1);
    
    for si1 = 1 :r
        % part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si1,r);
        
        Util= U(:,1:si1);
        Sigtil = S(1:si1,1:si1);
        Vtil = V(:,1:si1);
        
        all=[ Xd_p; Util'*X_p;out]*pinv([Xd; Sigtil*Vtil';inp]);
        
        approxA{si1} = all(1:size(Xd,1)+si1,1:size(Xd,1)+si1);
        approxB{si1} = all(1:size(Xd,1)+si1,size(Xd,1)+si1+1:end);
        approxC{si1}= all(size(Xd,1)+si1+1:end, 1:size(Xd,1)+si1);
        approxD{si1} = all(size(Xd,1)+si1+1:end, size(Xd,1)+si1+1:end);
        
        sys_red{si1} = ss(approxA{si1},approxB{si1},approxC{si1},approxD{si1},2);
        
        %same as before
        [FITje,fig1,x] = evaluatemodel(sys_red,si1,Inputs,Outputs,FITje,'identification',x,states,U,Deterministic,method);
        if ~isempty(fig1)
        warning off
        export_fig(fig1,strcat(dirdmdident,'/image',num2str(10000+si1)),'-nocrop','-m2')
        warning on
        close all
        end
    end
    
    [fig200]=VAFpermodes(FITje,r,{});
    warning off
    export_fig(fig200,strcat(dirdmdident,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
    warning on
    
    
elseif method == 4
    %% Professor Wingerden Least Square Solution for state space problem
    
    dirdmd = 'DMDresults_Wing';
    dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir')
        mkdir(dirdmd);
    end
    
    % Take singular value decomposition of X with rank r
    % X ~ USV*
    % U belongs to set C with size nxr
    % S belongs to set C with sixe rxr
    % V belongs to set C iwth size mxr, and * denotes the conjugate transpose
    % r is the rank of the reduced SVD approximation to X
    
    % left singular vectors U are POD modes
    % Columns of U are orthonormal, so U*U=I and V*V=I
    
    
    %Projection of full state space onto POD modes
    %Make use of SVD decompositoin in this phase
    
    % U*X' ~ U*A(USV*) + U*BU
    % Y  ~ C(USV*) + DU
    
    %A new state X^= U*X=U*USV*=SV*
    % ?  = U*A*U
    % B^ = U*B
    
    % U*X' ~ ?SV* + B^U
    
    % Matrixes ? and B^can now be found via least squares
    
    % || U*X' - [ ? B^][ SV* ] ||
    %                  [  U  ]
    
    % [ ? B^ ] = U*X [ SV* ]* x  [ SV*VS    SV*U_*  ] ^-1
    %                [ U_   ]    [ U_SV     U_U_*   ]
    
    %%%---%%%---%%%
    
    %including deterministic states, the matrix problem formulation will
    %be:
    
    % [ I 0  ] [ Xd ] ~ [ I 0  ] A [ Xd   ] + [ I 0  ] BU
    % [ 0 U* ] [ X' ] ~ [ 0 U* ]   [ USV* ]   [ 0 U* ]
    
    % with
    %
    % ? = [ I 0  ] A
    %     [ 0 U* ]
    
    % || [ I 0  ] [ Xd ]  -  ?[ Xd   ] - B^ U    ||
    % || [ 0 U* ] [ X' ]      [ USV* ]           ||
    
    % [ ? B^ ]
    
    %truncation/number of singular values for SVD decomposition
    
    Xd=[X1; X2; X3;X4];
    a=size(Xd);
    % states=[states];
    [Uo,So,Vo]=svds(X,r);
    U = blkdiag(eye(size(Xd,1)),Uo);
    QQ = [Xd; states];
    X = QQ(:,1:end-1);
    X_p=QQ(:,2:end);
    S=blkdiag(eye(size(Xd,1)),So);
    V=[Xd(:,1:end-1)' Vo];
    FITje=zeros(2, r+a(1));
    OMEGA={};
    %DAMPING={};
    
    Attt = cell((r+a(1)),1);
    At = cell((r+a(1)),1);
    Bt = cell((r+a(1)),1);
    Ctt = cell((r+a(1)),1);
    Ct = cell((r+a(1)),1);
    Dt = cell((r+a(1)),1);
    for si1=1:1:(r+a(1))
        
        
        Attt{si1}=U(:,1:si1)'*QQ(:,2:end)*[S(1:si1,1:si1)*V(:,1:si1)';...
            [U1(:,1:end-1)]]'/...
            ([S(1:si1,1:si1)*V(:,1:si1)'*V(:,1:si1)*S(1:si1,1:si1) S(1:si1,1:si1)*V(:,1:si1)'* ...
            [U1(:,1:end-1)]'; ...
            [U1(:,1:end-1)]*V(:,1:si1)*S(1:si1,1:si1) [U1(:,1:end-1)]*[U1(:,1:end-1)]' ] ); %#ok<USENS,NBRAK> U1 should be defined
        
        At{si1}=Attt{si1}(:,1:si1);
        Bt{si1}=Attt{si1}(:,si1+1:end);
        Ctt{si1}=[Y1(:,1:end-1);Y2(:,1:end-1)]*pinv(([S(1:si1,1:si1)*V(:,1:si1)'; U1(:,1:end-1)])); %#ok<USENS>
        
        Ct{si1} = Ctt{si1}(:,1:si1);
        Dt{si1} = Ctt{si1}(:,si1+1:end);
        
        sys_red{si1}=ss(At{si1},Bt{si1},Ct{si1}, Dt{si1},2);
        
        [FITje,OMEGA,fig1] = evaluatemodel(sys_red,si1,Inputs,Outputs,FITje,OMEGA,'identification');
        
        warning off
        export_fig(fig1,strcat(dirdmd,'/image',num2str(10000+si1)),'-nocrop','-m2')
        warning on
        close all
    end
    
    fig200  = VAFpermodes(FITje,r,{});
    warning off
    export_fig(fig200,strcat(dirdmd,'/image',num2str(10000+si1+1)),'-nocrop','-m2')
    warning on
    
end



end

