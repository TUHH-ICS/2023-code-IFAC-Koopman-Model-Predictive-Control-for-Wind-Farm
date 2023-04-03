function fig1 = plotEDMDinputsEffWindUpdate(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,nx,strVal,Pvec,displayOff,Vinf)
% plotEDMDinputsEffWind plots the simulated outputs of the eDMD estimations
% and prints them into the directory given in dirFig.

if nargin < 6
    Vin = 0;
end
if nargin < 7
    nx = 0;
end

if nargin < 8
    strVal = 'Val.';
end

if nargin < 9
    displayOff = 0;
end

if nargin < 10
    Vinf = [];
end

if nx > 0
    strK =  sprintf(', $n_K$ = %d', nx);
else
    strK = '';
end

fig1 = figure;
posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:2),posDefault(3),posDefault(4)*1.2]);
cl = colormap('lines');
lw = 1;
fs = 10;
kplot = 60: length(ysim_val);

noSub = size(Inputs_val,1) + 2;
subplot(noSub,1,1)
plot(Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('$C_T$ [-]','interpreter','latex','Fontsize',fs);

if displayOff
    title(sprintf('%s: VAF Ur1: %d, Ur2: %d',...
    strVal,round(FITje_val)),'interpreter','latex','Fontsize',fs);
else
title(sprintf('%s Set%s:  VAF Ur1: %d, Ur2: %d',...
    strVal,strK,round(FITje_val)),'interpreter','latex','Fontsize',fs);
end


strLeg = sprintf('WT%d,',1:2);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast');
legend('boxoff')
set(gca,'fontsize', fs)

if noSub > 3
    subplot(noSub,1,2)
    if isempty(Vinf)
    plot(Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    else
         plot(Vinf(kplot),'k','LineWidth',lw); hold on; plot(Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:));
    end
    grid on; axis tight
    if Vin && ~isempty(Vinf)
        ylabel('$V$ [m/s]','interpreter','latex','Fontsize',fs);
        legend('$V_{\inf}$','$V_1$','interpreter','latex','Fontsize',fs,'orientation', 'horizontal','Location','NorthWest');
        legend('boxoff')
    elseif Vin
        ylabel('Input $V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    else
        ylabel('Input $\Phi$ [deg]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
end

subplot(noSub,1,noSub-2)
plot(Deterministic_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(ysim_val(1,kplot),'Color',[0,0,0.6],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight;
pos = axis; axis([pos(1:3), pos(4)+0.3]);
%xlabel('k [-]','interpreter','latex')
ylabel('$U_{r1}$ [m/s]','interpreter','latex')
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
%  0.6155    0.4802    0.2702    0.0397
lgh.Position = [lgh.Position(1) + 0.02,lgh.Position(2:4)];
legend('boxoff')
set(gca,'fontsize', fs)

subplot(noSub,1,noSub-1)
plot(Deterministic_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)); %1
hold on;
plot(ysim_val(2,kplot),'Color',[0.6, 0, 0.0],'LineWidth',lw,'LineStyle','--') %3
grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
ylabel('$U_{r2}$ [m/s]','interpreter','latex','fontsize', fs)
lgh = legend({'Real','Est.'},'Location','southeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
lgh.Position = [0.2595    0.3585    0.2702    0.0397];
legend('boxoff')
set(gca,'fontsize', fs)

subplot(noSub,1,noSub)
plot(min((Pvec-1),1),'k','LineWidth',lw); %1

grid on; axis tight
pos = axis; axis([pos(1:3), pos(4)+0.3]);
xlabel('k [-]','interpreter','latex')
ylabel('$Update$ ','interpreter','latex')
axis tight;

%legend('P','Location','northeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
%legend('boxoff')
set(gca,'fontsize', fs)

%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
print(gcf,fullfile(dirFig,['ValP',strVal1]), '-dpng');
print(gcf,fullfile(dirFig,['ValP',strVal1]), '-depsc');