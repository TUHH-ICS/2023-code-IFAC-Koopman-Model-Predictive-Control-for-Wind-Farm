function fig1 = plotDMDoutputsFP(stateName,Outputs_val,ysim_val,FITje_val,FITje,Outputs,ysim,Deterministic_val,Deterministic)
%% plots for validation data
    fig1 = figure(1);
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
    fig1 = figure(2);
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
