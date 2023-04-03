function [sys_red,FITje,Xd,Xd_p,x,FITje_val,fig1,fig2,xo,Koop] = ...
    eDynamicmodedecomposition(states,statesvalidDMD,Inputs,Outputs,...
    Deterministic,Deterministic_val,Inputs_val,Outputs_val,method,...
    maindir,dirdmd,dt,stateName,dirFig,yawmode)

if nargin <14
    yawmode = 0;
end
%[sys_red,FITje,U,S,V,X,X_p,x]
%eDynamicmodedecomposition aims to build a reduced order model from the states,
% input/output information and deterministic states gathered in the
% simulation and resampled
nx = size(states,1);%for States Ur1,Ur2 and lifted states
nu = size(Inputs,1);
ny = size(Outputs,1);

Xd     = Deterministic(:,1:end-1);
Xd_p   = Deterministic(:,2:end);
NlUpdate = 1;

if method == -2
    if yawmode == 0
        
        X      = states(:,1:end-1);
        X_p    = states(:,2:end);
        out    = Outputs(:,1:end-1);
        inp    = Inputs(:,1:end-1);
        
        Psi = [X;inp];
        Koop.G = Psi* Psi';
        Koop.A = X_p* Psi';
        
        all1 = Koop.A*pinv(Koop.G); %X_p * pinv([X;inp]);%for States Ur1,Ur2 and lifted states
        approxA = all1(1:nx,1:nx);
        approxB = all1(1:nx,nx+1:end);
        approxC = [eye(2),zeros(2,nx-2)];
        approxD = zeros(2,nu);
        
        sys_red = ss(approxA,approxB,approxC,approxD,dt);
        %xo = Xd(:,1);%for States Ur1,Ur2
        xo = X(:,1);%for States Ur1,Ur2 and lifted states
        [ysim,~,xsim] = lsim(sys_red, Inputs',[],xo);
        x= xsim;
        FITje(:,1) = vaf(Deterministic',ysim);
        xo_val = statesvalidDMD(:,1)';%for States Ur1,Ur2 and lifted states
        [ysim_val,~,~] = lsim(sys_red, Inputs_val',[],xo_val);
        FITje_val = vaf(Deterministic_val',ysim_val(1:end,:));
        sys_red.StateName = strsplit(stateName,';');
        sys_red.OutputName = {'Ur1';'Ur2'};
        if max(size(sys_red)) == 3
            sys_red.InputName  =  {'Ct1';'Ct2';'phi1'};
        end
        %sys_red{1,
        
        %% plots for validation data
        fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val, FITje_val,dirFig);
    else
        X      = states(:,1:end-1);%states(:,1:end-1);
        X_p    = states(:,2:end);
        out    = Outputs(:,1:end-1);
        inp    = Inputs(:,1:end-1);
        
        Psi = [X;inp];
        Koop.G = Psi* Psi';
        Koop.A = [X_p;out]* Psi';
        
        all1 = Koop.A*pinv(Koop.G); %X_p * pinv([X;inp]);%for States Ur1,Ur2 and lifted states
        approxA = all1(1:nx,1:nx);
        approxB = all1(1:nx,nx+1:end);
        approxC = all1(nx+1:end,1:nx);
        approxD = all1(nx+1:end,nx+1:end);%[eye(2),zeros(2,nx-2)];
        
        sys_red = ss(approxA,approxB,approxC,approxD,dt);
        %xo = Xd(:,1)';%for States Ur1,Ur2
        xo = X(:,1);%for States Ur1,Ur2 and lifted states
        if NlUpdate == 1
            for i = 1:length(Inputs)
                xsim(:,i+1) =  approxA*xo(:,i)+approxB*Inputs(:,i);
                ysim(:,i) =  approxC*xo(:,i)+approxD*Inputs(:,i);
                xo(:,i+1)=[xsim(1,i+1),xsim(2,i+1),xsim(1,i+1)^2,xsim(2,i+1)^2,xsim(1,i+1)^3,xsim(2,i+1)^3]';
            end
            x= xo;
        else
            [ysim,~,xsim] = lsim(sys_red, Inputs',[],X(:,1));
            x= xsim;
        end
        FITje(:,1) = vaf(Outputs',ysim);
        xo_val = statesvalidDMD(:,1);%for States Ur1,Ur2 and lifted states
        if NlUpdate == 1
            for i = 1:length(Inputs_val)
                xsim_val(:,i+1) =  approxA*xo_val(:,i)+approxB*Inputs_val(:,i);
                ysim_val(:,i) =  approxC*xo_val(:,i)+approxD*Inputs_val(:,i);
                xo_val(:,i+1)=[xsim_val(1,i+1),xsim_val(2,i+1),xsim_val(1,i+1)^2,xsim_val(2,i+1)^2,xsim_val(1,i+1)^3,xsim_val(2,i+1)^3]';
            end
        else
            [ysim_val,~,~] = lsim(sys_red, Inputs_val',[],statesvalidDMD(:,1));
        end
        FITje_val = vaf(Outputs_val',ysim_val(:,1:end));
    end
    sys_red.StateName = strsplit(stateName,';');
    sys_red.OutputName = {'PT1';'PT2';'FT1';'FT2'};
    if max(size(sys_red)) == 3
        sys_red.InputName  =  {'Ct1';'Ct2';'phi1'};
    end
    %sys_red{1,
    
    %% plots for validation data
    if NlUpdate == 1
        fig2 = plotEDMDinputsEffWind(ysim',Inputs,Outputs, FITje,dirFig,yawmode,'Id');
        fig1 = plotEDMDinputsEffWind(ysim_val',Inputs_val,Outputs_val, FITje_val,dirFig,yawmode,'val');
    else
        fig2 = plotEDMDinputsEffWind(ysim,Inputs,Outputs, FITje,dirFig,yawmode,'Id');
        fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Outputs_val, FITje_val,dirFig,yawmode,'val');
    end
    
    %plotEDMDinputsEffWind(ysim,Inputs,Deterministic, FITje,dirFig)
    
    
%end



elseif method ==-1
    X      = states(:,1:end-1);
    X_p    = states(:,2:end);
    
    out    = Outputs(:,1:end-1);
    inp    = Inputs(:,1:end-1);
    
    Xd     = Deterministic(:,1:end-1);
    Xd_p   = Deterministic(:,2:end);
    dirdmd = 'eDMDresults_DMD';
    dirdmd = fullfile(maindir,dirdmd);
    if ~exist(dirdmd,'dir')
        mkdir(dirdmd);
    end
    %nx = size(Xd,1);%for States Ur1,Ur2
    nx = size(X,1);%for States Ur1,Ur2 and lifted states
    nu = size(inp,1);
    ny = size(out,1);
    %all1=[Xd_p;out]*pinv([Xd;inp]);%for States Ur1,Ur2
    all1=[X_p;out]*pinv([X;inp]);%for States Ur1,Ur2 and lifted states
    approxA = all1(1:nx,1:nx);
    approxB = all1(1:nx,nx+1:end);
    %approxC = all1(nx+1:end,1:nx); %[zeros(1,nx2), 1/nx2 * ones(1,nx2), zeros(1,si1)];
    approxC = [eye(2),zeros(2,nx-2);all1(nx+1:end,1:nx)];
    approxD = [zeros(2,nu);all1(nx+1:end,nx+1:end)];
    %approxD = all1(nx+1:end,nx+1:end);
    sys_red = ss(approxA,approxB,approxC,approxD,dt);
    %xo = Xd(:,1);%for States Ur1,Ur2
    xo = X(:,1);%for States Ur1,Ur2 and lifted states
    [ysim,~,xsim] = lsim(sys_red, Inputs',[],xo);
    x= xsim;
    %FITje(:,1) = vaf(Outputs',ysim);
    FITje(:,1) = vaf([Deterministic;Outputs]',ysim);
    %xo_val = [Deterministic_val(:,1)]';%for States Ur1,Ur2
    xo_val = statesvalidDMD(:,1)';%for States Ur1,Ur2 and lifted states
    [ysim_val,~,~] = lsim(sys_red, Inputs_val',[],xo_val);
    %FITje_val = vaf(Outputs_val',ysim_val(1:end,:));
    FITje_val = vaf([Deterministic_val;Outputs_val]',ysim_val(1:end,:));
    
    %% plots for validation data
    fig1=figure(1);
    fig1.Visible='on';
    fig1.WindowState = 'maximized';
    posDefault = [805.0000  272.2000  560.0000  420.0000];
    set(gcf, 'position', [posDefault(1:3),posDefault(4)*3]);
    %sgtitle(['U_1.^2;U_2.^2;U_1.^3;U_2.^3;\DeltaU_1;\DeltaU_2;\DeltaU_1.^2;\DeltaU_2.^2;\DeltaU_1.^3;\DeltaU_2.^3'])
    %set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','on');
    sgtitle(stateName)% for states U1 U2
    %sgtitle('States:[U_1,U_2,U^2_1,U^2_2,U^3_1, U^3_2,\DeltaU_1,\DeltaU_2,\DeltaU^3_1,\DeltaU^3_2,\DeltaU_1U_1,\DeltaU^2_2U2,\DeltaU^3_1U1,\DeltaU^3_2U2]')
    %,U^2_1, U^2_2,U^3_1, U^3_2,,\DeltaU^2_1,\DeltaU^2_2,\DeltaU_2U_2,\DeltaU^2_1U1,
    subplot(3,2,1)
    plot(Outputs_val(1,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim_val(:,3),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1 ]')
    ylabel(' Power [MW]')
    title({['P_{g1} VAF of ',num2str(round(FITje_val(3,1),2)),' %. '...
        'VAF on Ident. Data ',num2str(round(FITje(3,1),2)),' %.']})
    legend({'Real Generator Power','Model Output'},'Location','southwest')%'bestoutside','Orientation','horizontal')
    
    
    legend('boxoff')
    %set(gcf, 'position', [posDefault(1:3),posDefault(4)*10]);
    %set(gca,'fontsize', 12)
    
    subplot(3,2,2)
    plot(Outputs_val(2,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim_val(:,4),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1]')
    ylabel(' Power [MW]')
    title({['P_{g2} VAF of ',num2str(round(FITje_val(4,1),2)),' %. ']...
        ['VAF on Ident. Data ' ,num2str(round(FITje(4,1),2)),' %.']})
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    %set(gca,'fontsize', 12)
    subplot(3,2,3)
    plot(Outputs_val(3,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim_val(:,5),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1]')
    ylabel(' Thrust Force [kN]')
    title({['Ft_ VAF of ',num2str(round(FITje_val(5,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(5,1),2)),' %.']})
    legend({'Real Thrust Force','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    %set(gca,'fontsize', 12)
    subplot(3,2,4)
    plot(Outputs_val(4,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim_val(:,6),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1]')
    ylabel(' Thrust Force [kN]')
    title({['Ft_{2} VAF of ',num2str(round(FITje_val(6,1),2)),' %. ']...
        ['VAF on Ident. Data ' ,num2str(round(FITje(6,1),2)),' %.']})
    legend({'Real Thrust Force','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    %set(gca,'fontsize', 12)
    %saveas (gcf, 'dirdmd\statesUr.jpg');
    subplot(3,2,5)
    plot(Deterministic_val(1,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim_val(:,1),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1 ]')
    ylabel(' Wind velocity [m/s]')
    title({['U_ VAF of ',num2str(round(FITje_val(1,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(1,1),2)),' %.']})
    legend({'Real','Model Output'},'Location','bestoutside','Orientation','horizontal')
    
    
    legend('boxoff')
    %set(gca,'fontsize', 12)
    subplot(3,2,6)
    plot(Deterministic_val(2,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim_val(:,2),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1 ]')
    ylabel(' Wind velocity [m/s]')
    title({['U_{r2} VAF of ',num2str(round(FITje_val(2,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(2,1),2)),' %.']})
    legend({'Real','Model Output'},'Location','bestoutside','Orientation','horizontal')
    
    
    legend('boxoff')
    set(gca,'fontsize', 12)
    sys_red.StateName = strsplit(stateName,';');
    sys_red.OutputName = {'Ur1';'Ur2';'PT1';'PT2';'FT1';'FT2'};
    sys_red.InputName  =  {'Ct1';'Ct2'};
    %sys_red{1,1}.Xint = xo;
    
    %%
    % plots for Identification data
    fig1=figure(2);
    fig1.Visible='on';
    %sgtitle(['U_1.^2;U_2.^2;U_1.^3;U_2.^3;\DeltaU_1;\DeltaU_2;\DeltaU_1.^2;\DeltaU_2.^2;\DeltaU_1.^3;\DeltaU_2.^3'])
    set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','on');
    sgtitle(stateName)% for states U1 U2
    %sgtitle('States:[U_1,U_2,U^2_1,U^2_2,U^3_1, U^3_2,\DeltaU_1,\DeltaU_2,\DeltaU^3_1,\DeltaU^3_2,\DeltaU_1U_1,\DeltaU^2_2U2,\DeltaU^3_1U1,\DeltaU^3_2U2]')
    %,U^2_1, U^2_2,U^3_1, U^3_2,,\DeltaU^2_1,\DeltaU^2_2,\DeltaU_2U_2,\DeltaU^2_1U1,
    subplot(3,2,1)
    plot(Outputs(1,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim(:,3),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1 ]')
    ylabel(' Power [MW]')
    title({['P_{g1} VAF of ',num2str(round(FITje_val(3,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(3,1),2)),' %.']})
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    
    
    legend('boxoff')
    
    set(gca,'fontsize', 12)
    
    subplot(3,2,2)
    plot(Outputs(2,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim(:,4),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1]')
    ylabel(' Power [MW]')
    title({['P_{g2} VAF of ',num2str(round(FITje_val(4,1),2)),' %. ']...
        ['VAF on Ident. Data ' ,num2str(round(FITje(4,1),2)),' %.']})
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 12)
    subplot(3,2,3)
    plot(Outputs(3,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim(:,5),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1]')
    ylabel(' Thrust Force [kN]')
    title({['Ft_ VAF of ',num2str(round(FITje_val(5,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(5,1),2)),' %.']})
    legend({'Real Thrust Force','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 12)
    subplot(3,2,4)
    plot(Outputs(4,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim(:,6),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1]')
    ylabel(' Thrust Force [kN]')
    title({['Ft_{2} VAF of ',num2str(round(FITje_val(6,1),2)),' %. ']...
        ['VAF on Ident. Data ' ,num2str(round(FITje(6,1),2)),' %.']})
    legend({'Real Thrust Force','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 12)
    %saveas (gcf, 'dirdmd\statesUr.jpg');
    subplot(3,2,5)
    plot(Deterministic(1,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim(:,1),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1 ]')
    ylabel(' Wind velocity [m/s]')
    title({['U_ VAF of ',num2str(round(FITje_val(1,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(1,1),2)),' %.']})
    legend({'Real','Model Output'},'Location','bestoutside','Orientation','horizontal')
    
    
    legend('boxoff')
    set(gca,'fontsize', 12)
    subplot(3,2,6)
    plot(Deterministic(2,:),'b','LineWidth',1.6'); %1
    hold on;
    plot(ysim(:,2),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    axis tight
    xlabel('Time instant [dt = 1 ]')
    ylabel(' Wind velocity [m/s]')
    title({['U_{r2} VAF of ',num2str(round(FITje_val(2,1),2)),' %. ']...
        ['VAF on Ident. Data ',num2str(round(FITje(2,1),2)),' %.']})
    legend({'Real','Model Output'},'Location','bestoutside','Orientation','horizontal')
    
    else
        
        %define necessary matrices for DMD
        X      = states(:,1:end-1);
        X_p    = states(:,2:end);
        inp    = Inputs(:,1:end-1);
        
        %% eDMD Output is computed as the second part of the states
        dirdmd = fullfile(maindir,dirdmd);
        if ~exist(dirdmd,'dir')
            mkdir(dirdmd);
        end
        
        nx = size(X,1);
        nx2 = nx/2;
        nu = size(inp,1);
        
        Out1 = X;
        M1 = movmean(X(1,:),500);
        M2 = movmean(X(2,:),500);
        MP1 = movmean(X(1,:),500);
        MP2 = movmean(X(2,:),500);
        diff1 = X(1,:)-M1;
        diff2 = X(2,:)-M2;
        diff1_p = X_p(1,:)-MP1;
        diff2_p = X_p(2,:)-MP2;
        Xaug1 = [X.^2; X.^3;X(1,:).*X(2,:);diff1 ;diff2 ;diff1.^2;diff2.^2];
        X_p_aug1 = [X_p.^2; X_p.^3;X_p(1,:).*X_p(2,:);diff1_p ;diff2_p ;diff1_p.^2;diff2_p.^2];
        %nxAug = size(Xaug,1);
        
        [U,S,V] = svd(Xaug1);
        
        %Out1a = mean(X);
        
        r = size(U,1);
        FITje = NaN(r,2);
        
        approxA = cell(r,1);
        approxB = cell(r,1);
        approxC = cell(r,1);
        approxD = cell(r,1);
        sys_red = cell(r,1);
        x = cell(r,1);
        
        all1=X_p*pinv([X;inp]);
        approxA = all1(1:nx,1:nx);
        approxB = all1(1:nx,nx+1:end);
        approxC=[eye(nx)]; %[zeros(1,nx2), 1/nx2 * ones(1,nx2), zeros(1,si1)];
        approxD = zeros(nx,nu);
        sys_red = ss(approxA,approxB,approxC,approxD,dt);
        xo = states(:,1);
        [ysim,~,xsim] = lsim(sys_red, Inputs',[],xo);
        FITje(1,:) = vaf(Out1',ysim(1:end-1,:));
        x= xsim;
        for si1 = 1 : r
            % part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si1,r);
            
            Util= U(:,1:si1);
            Sigtil = S(1:si1,1:si1);
            Vtil = V(:,1:si1);
            
            
            % sys_red{si1} = ss(approxA{si1},approxB{si1},approxC{si1},approxD{si1},dt);
            X_p_aug = [X_p; Util'*X_p_aug1];
            Xaug = [X; Sigtil*Vtil'];
            all= X_p_aug * pinv([Xaug;inp]);
            %all=[X_p; Util'*X_p_aug1;Out1a]*pinv([X; Sigtil*Vtil';inp]);
            
            approxA{si1+1} = all(1:nx+si1,1:nx+si1);
            approxB{si1+1} = all(1:nx+si1,nx+si1+1:end);
            approxC{si1+1} =[eye(nx),zeros(nx,si1)]; %[zeros(1,nx2), 1/nx2 * ones(1,nx2), zeros(1,si1)];
            approxD{si1+1} = zeros(nx,nu);
            
            sys_red{si1+1} = ss(approxA{si1+1},approxB{si1+1},approxC{si1+1},approxD{si1+1},dt);
            
            
            xo = [states(:,1);  Util'*Xaug1(:,1);];
            [ysim,~,xsim] = lsim(sys_red{si1+1}, Inputs',[],xo);
            FITje(si1+1,:) = vaf(Out1',ysim(1:end-1,:));
            x{si1+1}= xsim;
            
        end
        
end
