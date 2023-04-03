close all; clear; clc;

%% Get current folder, data folder
try
    pwd1 = mfilename('fullpath');
catch
    pwd1 = pwd;
end

mainDir =  fileparts(fileparts(fileparts(fileparts(pwd1))));
filename = 'Vinf8dot0_sowfa_2turb_yaw_steps_solArray.mat';
mainDataDir = 'DataT2OLWFSim';
dataDir = 'Vinf8dot0_OL_Ct';

dirFigName = 'figPaper';
dirFig = fullfile(mainDir,dirFigName);

if ~isdir(dirFig) %#ok<ISDIR>
    mkdir(dirFig);
end

%% Load data
% filename = 'Vinf8dot0_sowfa_2turb_yaw_steps_solArray.mat';
% dataDir =  fullfile('DataT2OLWFSim','Vinf8dot0_OL_Ct_5Da');
% filename = 'Vinf8dot0_sowfa_2turb_yaw_steps_solArray.mat';
load(fullfile(mainDir,mainDataDir,dataDir,filename));

filenamepng = 'WFSimAnimation';
kPlotVec = (500:500:length(sol_array))-1;
angleVec = 0:5:35;
for idx = 1: length(kPlotVec)
    aK = kPlotVec(idx);
    sol = sol_array(aK);
    hfig(idx) = figure(aK); %#ok<SAGROW>
    pos = get(hfig(idx),'Position');
    set(hfig(idx),'Position',[pos(1:3),pos(4)*1.2]);
    WFSim_animation(Wp,sol,hfig(idx));
    
    filenamepng1 = sprintf('%s_Angle%02d',filenamepng,angleVec(idx));
    % saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
    print(gcf,fullfile(dirFig,filenamepng1), '-dpng');
    print(gcf,fullfile(dirFig,filenamepng1), '-depsc');
    
end


% hfig = figure; %#ok<SAGROW>
% pos = get(hfig,'Position');
% set(hfig,'Position',[pos(1:3),pos(4)*1.2]);
% 
% for idx = 1700 : 1900 %length(sol_array)
%     %aK = kPlotVec(idx);
%     sol = sol_array(idx);
%     WFSim_animation(Wp,sol,hfig);
%     F(idx) = getframe(gcf);
%     
% end


dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1);
fs = 12;

%% Define parameters
%a1Vec = repmat(1/3,1, length(gammaVec))*0.9; %
ctVec = [0,15,30,45;0.79, 0.76,0.69,0.58]; %temp = interp1([0,15,30,45],[0.79, 0.76,0.69,0.58],0:45);
cpVec = [0,15,30,45; 0.5, 0.48,0.42,0.26]; %temp = interp1([0,15,30,45],[0.5, 0.48,0.4,0.26],0:45);

gammaVec = (0:5:35); % change
alphaStar =  2.32; %32; % 4*alpha = 4* 0.58 = 2.32 2.8
betaStar = 0.154; % 2*beta = 0.154 from jetstreams
b   = 0.022; % 0.075; % model parameters wake interaction 0.022
I0 = 0.05; %Wp.site.v_Inf/Wp.site.u_Inf; % Intensity: 0 for our case


rho = Wp.site.Rho; % air density
D   = Wp.turbine.Drotor;  % rotor diameter
r   = D/2; %rotor radius
A   = r^2 *pi;% rotor area

Vinf  =  Wp.site.u_Inf; % wind unaffected by wind farm
xDist = diff(Wp.turbine.Crx); % 5*D; % Distance between turbines
%c12 = ( D/ (D + 2* b * xDist ) )^2;

a2 = 1/3*0.95; % Induction factor WT2 for maximal power  1/9; %
Ct2 = 4*a2 *(1- a2); % eq. 6.2 thrust coefficient for max. power
a1Vec =linspace(1/3, 1/3.5, length(gammaVec))*0.95; %

A12Vec = nan(length(gammaVec),1);
A12ToA1 = nan(length(gammaVec),1);
Ct1Vec = nan(length(gammaVec),1);
P = nan(length(gammaVec),2);

% For centerline and wake expansion
tempXDist = xDist; %this is the currently used distance
nD = ceil(10*D);
xDistVec = 0 : nD;
deltaVec = nan(nD,length(gammaVec));
sigmaYVec = nan(nD,length(gammaVec));
sigmaZVec = D/sqrt(8)*ones(nD,length(gammaVec));
x0Vec = nan(length(gammaVec),1);

% b = 0.022; % model parameters wake expansion 0.022 (from paper), 0.075
% commonly used
dUvecCell = cell(length(gammaVec),1);

% Position turbine
nTx = floor(Wp.turbine.Crx(1)/Wp.mesh.Lx * Wp.mesh.Nx);
nTy = round(Wp.turbine.Cry(1)/Wp.mesh.Ly * Wp.mesh.Ny);


% Values open loop: wind grid u,v and pressure grid 
u = cell2mat(arrayfun(@(x)(x.u), sol_array,'UniformOutput', false));
v = cell2mat(arrayfun(@(x)(x.v), sol_array,'UniformOutput', false));

kk = length(u)/Wp.mesh.Ny;

utemp = reshape(u,Wp.mesh.Nx,Wp.mesh.Ny,kk);
utempT = squeeze(utemp(nTx,nTy,:));
vtemp = reshape(v,Wp.mesh.Nx,Wp.mesh.Ny,kk);
vtempT = squeeze(vtemp(nTx,nTy,:));

figure; plot(diff(vtempT))


for idx = 1 : length(gammaVec)
    gammaDeg = gammaVec(idx); % current yaw angle WT1
    a1 = a1Vec(idx);
    gamma =  gammaDeg/180*pi;
    
    % Calculate variables based on paper
    % eq.6.7 Thrust coefficient (approximation of eq. 6.1
    % Ct1 = 4*a1 *(1- a1*cos(gamma)); %* interp1(ctVec(1,:),ctVec(2,:),gammaDeg);
    Ct1 = 4*a1 * sqrt(1- a1*(2*cos(gamma)-a1)); %* interp1(ctVec(1,:),ctVec(2,:),gammaDeg)/ctVec(2,1);
    % eq. 6.12: Wake angle at rotor 30 -> 6 deg
    ThetaCo = - 0.3*gamma/cos(gamma) *(1- sqrt(1-Ct1*cos(gamma)));
    
    % Question phi = - sin(gamma)*Ct1/(2*(1+ b * xDist/D)^2): 14 deg >> 6 deg
    I = I0 +   0.0150 *(idx-1); %%abs(sol_array(kPlotVec(idx)).v(nTx,nTy)/sol_array(kPlotVec(idx)).u(nTx,nTy));
    Ivec(idx) = I; 
    % end of near wake area x0
    sqrt1minCt = sqrt(1-Ct1);
    x0 = cos(gamma)* (1 + sqrt1minCt) / ...
        (sqrt(2)*(alphaStar*I + betaStar*(1 - sqrt1minCt))) * D; % eq. 6.16/ 7.3
    x0Vec(idx) = x0;
    
    % Calculate centerline and wake expansion
    for idx2 = 1: length(xDistVec)
        xDist =  xDistVec(idx2);
        % Wake width
        if xDist <= x0
            deltaVec(idx2,idx) = ThetaCo * xDist;
            sigmaYVec(idx2,idx) = D * cos(gamma)/sqrt(8);
        else
            
            % wake defelection delta
            term1 = ThetaCo * x0;
            % 7.2 eq. for sigmaY and sigmaZ
            sigmaY = b * (xDist - x0) + D * cos(gamma)/sqrt(8);
            sigmaZ = b * (xDist - x0) + D/sqrt(8);
            
            sqrtCt = sqrt(Ct1);
            tempTerm = 1.6 * sqrt((8*sigmaY*sigmaZ)/(D^2*cos(gamma)));
            
            term2 = ThetaCo * D/14.7 * sqrt(cos(gamma)/(b^2*Ct1)) *...
                (2.9 + 1.3 *sqrtCt - Ct1) * ...
                log( (1.6 + sqrtCt)*(tempTerm - sqrtCt)/...
                ((1.6 - sqrtCt)*(tempTerm + sqrtCt)));
            
            deltaVec(idx2,idx) = term1 + term2;
            sigmaYVec(idx2,idx) = sigmaY;
            sigmaZVec(idx2,idx) = sigmaZ;
        end
    end
    
    % Get wake deflection and expansion at second turbine
    xDist = tempXDist;
    [~,idxD] = min(abs(xDistVec - xDist));
    d = deltaVec(idxD,idx);
    
    sigmaY  = sigmaYVec(idxD,idx); % wake widths at xDist;
    R = sigmaY;% wake widths at xDist;
    sigmaZ = sigmaZVec(idxD,idx);
    
     %Eq. 7.1 far wake turbine deficit: sigma
    fracCT = (Ct1*cos(gamma)) / (8*(sigmaY*sigmaZ/D^2));
    
    yvec = -r : r ; %d - R : d + R;
    zvec = -r:r;
    dUVec = nan(length(yvec),length(zvec));
    for idxZ = 1: length(zvec)
        zDist = zvec(idxZ);
        for idxY = 1: length(yvec)
            yDist = yvec(idxY);
            if (zDist^2 + yDist^2 < r^2)
                dUVec(idxY,idxZ) = (1- sqrt(1 - fracCT ))*exp(-0.5*( ((yDist - d)/sigmaY)^2 + (zDist/sigmaZ)^2));
            end
        end
    end
   
    dUvecCell{idx} = dUVec;
    
    dURWake = mean(dUVec(:),'omitnan');
    dWT2 = (1 - dURWake);
    
    deltaV = [1, dWT2]; %];
    % deltaV =  [0,(c12* a1 * A12ToA1Vec)];
    U = Vinf .* deltaV;   %(1 - 2 * deltaV).* %[cos(gamma), 1]; %disp(V) ToDo: This is only for two turbines
    Cp1 = 1; %interp1(cpVec(1,:),cpVec(2,:),gammaDeg)/cpVec(2,1);
    Cp = 4*[a1,a2].*([sqrt(1- a1*(2*cos(gamma)-a1)),(1- a2)]).^2;
    Cp = [Cp(1)*Cp1, Cp(2)]; %*0.6
    
    P(idx,:) = 1/2* rho * A* Cp .* (U .* [cos(gamma), 1]).^3;%  - [0, 5.4798e+05];
    Vvec(idx,:) = U;
    ThetaCoVec(idx) = ThetaCo;
end

% Plot 2D wake on second turbine
figure; 
for idx = 1:4
    subplot(2,2,idx); 
    pcolor(fliplr(dUvecCell{2*idx-1}')); shading flat; 
    cb = colorbar;  set(cb,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight; 
    set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
    title(sprintf('$\\gamma$ = %d deg',gammaVec(2*idx-1)))
end

% Plot Centerline
figure;
plot(xDistVec/D, deltaVec(:,1:7)/D,xDistVec/D, deltaVec(:,8)/D,'k',...
    'linewidth', 1);
axis tight;  grid on;
temp = regexp(sprintf('%d deg,', gammaVec),',','split');
legend(temp(1:end-1),'location','southwest','Interpreter','Latex','FontSize',fs);

xlabel('Normalized Distance from Turbine $x/d$ [-]');
ylabel('Normalized Wake Deflection $\delta/d$ [-]')
set(0,'DefaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','Latex','FontSize',fs)

strName = sprintf('CenterlineDeflection');
set(gcf,'Name', strName);
print(gcf,fullfile(dirFig,strName), '-dpng');
print(gcf,fullfile(dirFig,strName), '-depsc');


% Plot Centerline and wake centers
figure;
plot(xDistVec, deltaVec(:,1),'b', xDistVec, deltaVec(:,1) + sigmaYVec(:,1), 'b--',...
    xDistVec, deltaVec(:,1) - sigmaYVec(:,1), 'b-.');
hold on;
plot(xDistVec, deltaVec(:,5),'r', xDistVec, deltaVec(:,5) + sigmaYVec(:,5), 'r--',...
    xDistVec, deltaVec(:,5) - sigmaYVec(:,5), 'r-.');
axis tight;  grid on;


figure;
vec = [1,5];
lc = lines;
%lc = [unique(lc,'rows');zeros(1,3)];
for idx = vec
    plot(xDistVec, deltaVec(:,idx), xDistVec, deltaVec(:,idx) + sigmaYVec(:,idx), '--',...
        xDistVec, deltaVec(:,idx) - sigmaYVec(:,idx), '-.','Color',lc(idx,:));
    hold on;
end
axis tight;  grid on;
temp = regexp(sprintf('$\\delta$ %d deg, $\\delta$+$\\sigma$ %d deg, $\\delta$-$\\sigma$ %d deg,', repmat(gammaVec(vec),3,1)),',','split');
legend(temp(1:end-1),'location','southwest','Interpreter','Latex');

% Comparison WFSim/ Gaussian Wake model
Psum = sum(P,2);
[Pmax,idxPmax]= max(Psum);
Phi = cell2mat(arrayfun(@(x)(x.turbine.Phi),sol_array,'UniformOutput', false));
phi1 = Phi(1,:); phi2 = Phi(2,:);
Power = cell2mat(arrayfun(@(x)(x.turbine.power),sol_array,'UniformOutput', false));
PT1 = Power(1,:); PT2 = Power(2,:);
U = cell2mat(arrayfun(@(x)(x.turbine.Ur),sol_array,'UniformOutput', false));
Ur1 = U(1,:);  Ur2 = U(2,:);

Ur1Plot = Ur1(kPlotVec);
Ur2Plot = Ur2(kPlotVec);

PT1Plot = PT1(kPlotVec);
PT2Plot = PT2(kPlotVec);

phi1K = phi1(kPlotVec);
Pk = PT1(kPlotVec)+ PT2(kPlotVec);
vecIdx = (max(Pk) == Pk);

figure;
plot(phi1(kPlotVec),PT1(kPlotVec)/10^6,'bo-', phi1(kPlotVec),PT2(kPlotVec)/10^6,'ro-', ...
    phi1(kPlotVec),Pk/10^6,'ko-',phi1K(vecIdx),max(Pk)/10^6,'go','Linewidth',1); ...
    %legend('P WT1','P WT2', 'P T','max(P T)','Location','eastoutside');
axis tight; grid on; hold on;
set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
xlabel('yaw $\gamma$ [deg]','Interpreter','Latex','FontSize',fs);
ylabel('Power [MW]','Interpreter','Latex','FontSize',fs);

legend('$P1_{WFSim}$','$P2_{WFSim}$','$PT_{WFSim}$','max($PT_{WFSim}$)',...
    'Interpreter','Latex','Location','EastOutside');

plot(gammaVec,P(:,1)/10^6,'bs--',gammaVec,P(:,2)/10^6,'rs--');
title(sprintf('WT1-WT2 [m]: %d, V [m/s]: %d', round(xDist), round(Vinf)));
plot(gammaVec,Psum/10^6,'ks--',gammaVec(idxPmax),Pmax/10^6,'gs')

strName = sprintf('PowerVsYawWFSim');
set(gcf,'Name', strName);
print(gcf,fullfile(dirFig,strName), '-dpng');
print(gcf,fullfile(dirFig,strName), '-depsc');


hl = legend({'$P1_{WFSim}$','$P2_{WFSim}$','$PT_{WFSim}$','max($PT_{WFSim}$)',...
    '$P1_{Gauss}$','$P2_{Gauss}$','$PT_{Gauss}$','max($PT_{Gauss}$)'}',...
    'Interpreter','Latex','Location','EastOutside');

%sprintf('max P: $\\gamma$ = %d [deg]', phi1K(vecIdx)),...
strName = sprintf('PowerVsYaw');
set(gcf,'Name', strName);
print(gcf,fullfile(dirFig,strName), '-dpng');
print(gcf,fullfile(dirFig,strName), '-depsc');


%% Run image sweep again but pass estimated centerline
%cosGamma = cos(angleVec./180*pi);


for idx = 1: length(kPlotVec)
    aK = kPlotVec(idx);
    sol = sol_array(aK);
    hfig(idx) = figure(aK+1); %#ok<SAGROW>
    pos = get(hfig(idx),'Position');
    set(hfig(idx),'Position',[pos(1:3),pos(4)*1.2]);
    
    idxPlot = xDistVec <= diff(Wp.turbine.Crx);
    idxPlotSigmaNaN = xDistVec(idxPlot) < (x0Vec(idx)-5);
    delta = deltaVec(idxPlot,idx);
    sigmaY = sigmaYVec(idxPlot,idx);
    sigmaY(idxPlotSigmaNaN) = NaN;
    xDistPlot = xDistVec(idxPlot);
    
    WFSim_animation(Wp,sol,hfig(idx), delta,sigmaY, xDistPlot);
    
    filenamepng1 = sprintf('%s_AngleCenter%02d',filenamepng,angleVec(idx));
    % saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
    print(gcf,fullfile(dirFig,filenamepng1), '-dpng');
    print(gcf,fullfile(dirFig,filenamepng1), '-depsc');
    
end


set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);

return; 


%% test plot centerline
Dr = Wp.turbine.Drotor; %#ok<UNRCH>
Cry = Wp.turbine.Cry;
Crx = Wp.turbine.Crx;

input.phi = 35;

turb_coord = .5*Dr*exp(1i*input.phi*pi/180);
kk = 1;
Qy = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
mIdxQ = round(length(Qy)/2);
Qy1 = Qy(mIdxQ);
Qx = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
Qx1 = Qx(mIdxQ);
%
% figure; plot(-deltaVec(:,8) + Qy1, xDistVec + Qx1)
%
ldxx   = Wp.mesh.ldxx;
ldyy   = Wp.mesh.ldyy;
ldxx2  = Wp.mesh.ldxx2;
yline  = Wp.mesh.yline;

u_Inf  = Wp.site.u_Inf;

N      = Wp.turbine.N;
Cry    = Wp.turbine.Cry;
Crx    = Wp.turbine.Crx;
input  = sol.turbInput;

idxX = 5:size(sol.u,1)-4;
idxY = 5:size(sol.u,2)-4;

%     %% Plot u velocity flow component
figure;
contourf(ldyy(1,idxY ),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.2),'Linecolor','none');

idxX = 5:size(sol.u,1)-4;
idxY = 5:size(sol.u,2)-4;

%% Plot u velocity flow component
contourf(ldyy(1,idxY ),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.2),'Linecolor','none');
colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.04]);  hold on; hc = colorbar;
set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
axis equal; axis tight;
% plot(-deltaVec(:,8) + Qy1, xDistVec + Qx1)
plot(-deltaVec(:,8)+ Qy1 ,xDistVec+ Qx1, 'g');
plot(-deltaVec(:,8)+ Qy1 + sigmaYVec(:,8),xDistVec+ Qx1, 'g--');
plot(-deltaVec(:,8)+ Qy1 - sigmaYVec(:,8),xDistVec+ Qx1, 'g--');

set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);
