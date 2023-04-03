function [Wp] = layoutSet_sowfa_2turb_alm_turbl_Vinf()

Wp = struct('description','2 NREL 5MW turbines case, turbulent inflow, with NTM wind');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747],... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.40 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',8,... % Initial long. wind speed in m/s 8.0446
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.2,... % Mixing length in x-direction (m)
    'd_lower',73.3,... % Turbulence model gridding property
    'd_upper',601.9,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );



Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );


U = Wp.site.u_Inf;
Iref = 0.1400;
Lambda1 = 63; %63.0000;
T =  600000; % This was changed
t = 0 : 1:T;
%Vref = 50;

%% From Wind.m, ll 678-709

% NTM values
sigma1 = Iref*(0.75*U+5.6);
sigma2 = 0.8*sigma1;
Lu = 8.1*Lambda1;
Lv = 2.7*Lambda1;

% Wave number vector
L = U * T;
N = length(t);
m = ifftshift(-N/2:N/2-1);
k = 2*pi*m/L;

% Spectrum
Fu = sigma1^2 * 4*Lu/U./(1+6/(2*pi)*abs(k)*Lu).^(5/3);
Fv = sigma2^2 * 4*Lv/U./(1+6/(2*pi)*abs(k)*Lv).^(5/3);
n = randn([2,length(k)]) + sqrt(-1)*randn([2,length(k)]);
dZ = sqrt(2*pi*[Fu;Fv]/L) .* n;

% IFFT
u = N*real(ifft(dZ(1,:)));
v = N*real(ifft(dZ(2,:)));
u = (u -mean(u)) * sigma1/std(u) + U;
v = (v -mean(v)) * sigma2/std(v);

Wp.site.u_Inf = u;
Wp.site.v_Inf = v;





% Tuning notes '2turb_alm_turb' (Sep 6th, 2017):
% Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
% Retuned the turbulence model by hand on November 5th, 2018
end

