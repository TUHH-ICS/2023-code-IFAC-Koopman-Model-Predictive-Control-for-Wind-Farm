function turbInputSet = controlSet_sowfa_2turb_yaw_alm_combined2(Wp)
% controlSet_sowfa_2turb_yaw_steps combines controlSet_sowfa_2turb_yaw_alm_turbl
% and controlSet_sowfa_2turb_yaw_alm_uniform with added filtered white noise
% Testsignal for MIMO Koopman controller

if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
end

turbInputSet1 = controlSet_sowfa_2turb_yaw_alm_turbl;
turbInputSet2 = controlSet_sowfa_2turb_yaw_alm_uniform;

k1 = 1000;%875; %te = turbInputSet1.t(k1);
k2 = 500;%60; %t0 = turbInputSet1.t(k2);

 %turbInputSet.t = [turbInputSet1.t(1:k1), turbInputSet2.t(k2:end-1) +(te +1 -t0)]
    
% tempPhi(1,:) = [turbInputSet1.phi(1,1:k1), turbInputSet2.phi(1,k2:end)...
 %   turbInputSet1.phi(1,1:k1)+20, turbInputSet2.phi(1,k2:end)+20];
%tempPhi(2,:) = zeros(1,length(tempPhi(1,:)));
%tempPhi(2,:) = [turbInputSet1.phi(2,1:k1), turbInputSet2.phi(2,k2:end)...
   % turbInputSet1.phi(2,1:k1), turbInputSet2.phi(2,k2:end)];
    
%turbInputSet.t = 0: length(tempPhi)-1;

temp = [turbInputSet1.CT_prime(:,1:k1), ...
    turbInputSet2.CT_prime(:,k2:end),turbInputSet1.CT_prime(:,1:k1), ...
    turbInputSet2.CT_prime(:,k2:end),turbInputSet1.CT_prime(:,1:k1), ...
    turbInputSet2.CT_prime(:,k2:end),turbInputSet1.CT_prime(:,1:k1), ...
    turbInputSet2.CT_prime(:,k2:end)];
%l = round(length(temp)/12);
tempPhi(1,:) = [turbInputSet1.phi(1,1:k1), turbInputSet2.phi(1,k2:end)+10,...
    turbInputSet1.phi(1,1:k1)+25, turbInputSet2.phi(1,k2:end)+20,...
    turbInputSet1.phi(1,1:k1)+5, turbInputSet2.phi(1,k2:end)+15,...
    turbInputSet1.phi(1,1:k1)+10, turbInputSet2.phi(1,k2:end)];

phiRLim = 1;

for idx = 2:length(tempPhi(1,:)) % limit the yaw rate
    aDiff = tempPhi(1,idx) - tempPhi(1,idx-1);
    tempPhi(1,idx) = tempPhi(1,idx-1) + min(max(aDiff,-phiRLim),phiRLim);
end

%        30*ones(l,1);40*ones(l,1); 50*ones(l,1);60*ones(l,1); -10*ones(l,1);-20*ones(l,1);...
%        -30*ones(l,1);-40*ones(l,1); -50*ones(length(temp)-round(11*l),1)];
tempPhi(2,:) = zeros(1,length(tempPhi(1,:)));
turbInputSet.t = 0: length(tempPhi)-1;

tfD = c2d(tf(1,[2000,1]),1);
B = tfD.Numerator{:};
A = tfD.Denominator{:};

% X = randn(1,k1);
% X1 = randn(1,k1);
% Xfilt(1,:) = filter(B,A,X);
% Xfilt(2,:) = filter(B,A,X1);
% rnd1 = [Xfilt, zeros(2,length(turbInputSet2.CT_prime(:,k2:end)))];

%X = randn(1,k1+length(turbInputSet2.CT_prime(:,k2:end)))*2;
%X1 = randn(1,k1+length(turbInputSet2.CT_prime(:,k2:end)))*2;
X = randn(1,length(temp))*1;
X1 = randn(1,length(temp))*1;
Xfilt(1,:) = filter(B,A,X);
Xfilt(2,:) = filter(B,A,X1);
rnd1 = Xfilt;

tfD1 = c2d(tf(1,[30,1]),1);
B1 = tfD1.Numerator{:};
A1 = tfD1.Denominator{:};

X = (randn(1,length(temp)))*5;
X1 = (randn(1,length(temp)))*5;
Xfilt1(1,:) = filter(B1,A1,X);


% stdT = std(temp);
tmpCt  = min(max((temp + rnd1),0.2),2) ;
CTRLim = 0.01;
for idx = 2:length(tmpCt(:,:)) % limit the yaw rate
    aDiff = tmpCt(:,idx) - tmpCt(:,idx-1);
    tmpCt(:,idx) = tmpCt(:,idx-1) + min(max(aDiff,-CTRLim),CTRLim);
end


turbInputSet.CT_prime = tmpCt; 

tempPhi(1,:)  = tempPhi(1,:) + Xfilt1(1,:); 
for idx = 2:length(tempPhi(1,:)) % limit the yaw rate
    aDiff = tempPhi(1,idx) - tempPhi(1,idx-1);
    tempPhi(1,idx) = tempPhi(1,idx-1) + min(max(aDiff,-phiRLim),phiRLim);
end

turbInputSet.phi =  tempPhi;
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time

% For visualization:
% turbInputSet = controlSet_sowfa_2turb_yaw_steps;
% figure; subplot(2,1,1); plot(turbInputSet.CT_prime'); axis tight; grid on;
% subplot(2,1,2); plot(turbInputSet.phi'); axis tight; grid on;

if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
    error('Number of turbines in layout does not match your controlSet.');
end
end