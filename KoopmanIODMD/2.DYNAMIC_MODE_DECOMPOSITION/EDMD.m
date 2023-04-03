function [sys_red,FITje,U,S,V,method,X,X_p,Xd,diredmd,x] = EDMD(states, Inputs, Outputs, Deterministic,method,r,maindir)
% x(k+1) ~ A x(k)
% X' ~ AX

% X = [   |  |        |  ]
%     [   x1 x2 ... xm-1 ]
%     [   |  |        |  ]

% X'= [   |  |        |  ]
%     [   x2 x3 ...  xm  ]
%     [   |  |        |  ]

%define necessary matrices for DMD
diredmd = 'EDMDresults_DMD';
diredmd = fullfile(maindir,diredmd);
if ~exist(diredmd,'dir')
    mkdir(diredmd);
end
X      = states(:,1:end-1);
X_p    = states(:,2:end);

out      = Outputs(:,1:end-1);
inp      =Inputs(:,1:end-1);

Xd     = Deterministic(:,1:end-1);
Xd_p   = Deterministic(:,2:end);


plotOn = 1;


[FITje,fig1,x] = evaluatemodel(sys_red,si1,Inputs,Outputs,FITje,'identification',x,states,U,Deterministic,method,plotOn);

if ~isempty(fig1)
    warning off
    export_fig(fig1,strcat(dirdmdident,'/image',num2str(10000+si1)),'-nocrop','-m2')
    warning on
    close all
end
[fig200]=VAFpermodes(FITje,r,{});
warning off
export_fig(fig200,strcat(dirdmdident,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
warning on
