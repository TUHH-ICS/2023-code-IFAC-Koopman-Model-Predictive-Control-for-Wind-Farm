function fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val,FITje_val,dirFig,Vin,nx,strVal,displayOff,Vinf,strOut)
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

if nargin < 11
    strOut = 'U';
end

if nx > 0
    strK =  sprintf(', $n_K$ = %d', nx);
else
    strK = '';
end


if size(ysim_val,2) < size(ysim_val,1)
    ysim_val = ysim_val';
end


fig1 = figure;
cl = lines;
lw = 1;
fs = 10;
kplot = 60: length(ysim_val);

noSubIn = size(Inputs_val,1);
noSubOut =  min(size(ysim_val));
noSub = noSubOut + noSubIn -1;
subplot(noSub,1,1)
plot(Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
hold on;
plot(Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
grid on; axis tight; ylim([0 2+0.2])
ylabel('$C_T$ [-]','interpreter','latex','Fontsize',fs);

if displayOff && length(FITje_val) >1
    title(sprintf('%s:  VAF %sr1: %d, r2: %d',...
        strVal,strOut,round(FITje_val(1:2))),'interpreter','latex','Fontsize',fs);
elseif length(FITje_val)>1
    title(sprintf('%s Set%s:  VAF %sr1: %d, r2: %d',...
        strVal,strK,strOut,round(FITje_val(1:2))),'interpreter','latex','Fontsize',fs);
elseif displayOff
    title(sprintf('%s:  VAF %sT: %d',...
        strVal,strOut,round(FITje_val)),'interpreter','latex','Fontsize',fs);
else
    title(sprintf('%s Set%s:  VAF %sT: %d',...
        strVal,strK,strOut,round(FITje_val)),'interpreter','latex','Fontsize',fs);
end


strLeg = sprintf('WT%d,',1:2);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast');
%legend('boxoff')
set(gca,'fontsize', fs)

if noSubIn > 2
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
        %legend('boxoff')
    elseif (Vin && noSub == 4) || (~Vin && noSub == 3) 
         ylabel('$\gamma$ [deg]','interpreter','latex','Fontsize',fs);       
    else
        ylabel('$V_1$ [m/s]','interpreter','latex','Fontsize',fs);
    end
    set(gca,'fontsize', fs)
    
    if noSubIn > 3 %&& noSub ~= 4
        subplot(noSub,1,3)
        plot(Inputs_val(4,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
        ylabel('$V_1$ [m/s]','interpreter','latex','Fontsize',fs);
        grid on; axis tight
    end
    
end

if noSubOut == 1
    strAdd = {'_T$ [MW]'};
else
    strAdd = {'_{r1}$ [m/s]','_{r2}$ [m/s]'};
end

for idx = 1: noSubOut
    idxT = 2-mod(idx,2);
    subplot(noSub,1,noSub - (noSubOut-idx))
    plot(Deterministic_val(idx,kplot),'LineWidth',lw,'Color', cl(idxT,:)); %1
    hold on;
    plot(ysim_val(idx,kplot),'Color',[0,0,0],'LineWidth',lw,'LineStyle','--') %3
    grid on; axis tight;
    pos = axis; axis([pos(1:3), pos(4)+0.3]);

    %xlabel('k [-]','interpreter','latex')
    ylabel(['$',strOut,strAdd{idxT}],'interpreter','latex')
    legend({'Real','Est.'},'Location','northwest','interpreter','latex','fontsize',fs,'orientation','horizontal')
    %legend('boxoff')
    set(gca,'fontsize', fs)
    
end


%sys_red{1,1}.Xint = xo;
strVal1 = matlab.lang.makeValidName(strVal);
print(gcf,fullfile(dirFig,['Val',strVal1]), '-dpng');
print(gcf,fullfile(dirFig,['Val',strVal1]), '-depsc');