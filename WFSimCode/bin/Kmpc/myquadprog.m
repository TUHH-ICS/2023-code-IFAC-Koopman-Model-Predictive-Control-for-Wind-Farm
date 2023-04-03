function [X,fval] = myquadprog(H,f,A,B,varargin)
% myquadprog Quadratic programming based on QUADPROG with flexibility removed
%for speed-up
%
% Private function [X, fval] = ...
%    ipqpdense(full(H), f, A, B, Aeq, Beq, lb, ub, X0, flags, ...
%    options, defaultopt);
% is called 
%   X = QUADPROG(H,f,A,b) attempts to solve the quadratic programming
%   problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x <= b
%             x
%
%   X = QUADPROG(H,f,A,b,Aeq,beq) solves the problem above while
%   additionally satisfying the equality constraints Aeq*x = beq. (Set A=[]
%   and B=[] if no inequalities exist.)

defaultopt = struct( ...
    'Algorithm','interior-point-convex', ...
    'Diagnostics','off', ...
    'Display','final', ...
    'HessMult',[], ...
    'MaxIter',[], ...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
    'PrecondBandWidth',0, ...
    'ProblemdefOptions', struct, ...
    'TolCon',1e-8, ...
    'TolFun',[], ...
    'TolFunValue', [], ...
    'TolPCG',0.1, ...
    'TolX',100*eps, ...
    'TypicalX','ones(numberOfVariables,1)', ...
    'LinearSolver', 'auto', ...
    'ObjectiveLimit', -1e20 ...
    );

% Handle missing arguments
% Options are all default
options = defaultopt;
optimgetFlag = 'optimoptions';
computeFirstOrderOpt = false;

% Options setup
display = optimget(options,'Display',defaultopt,optimgetFlag);
detailedExitMsg = contains(display,'detailed');
verbosity = 0;
numberOfVariables = size(A,2);

Aeq = zeros(0,numberOfVariables);
Beq = zeros(0,1);
f = f(:);
B = B(:);
Beq = Beq(:);

%[X0,lb,ub] = checkbounds(X0,lb,ub,numberOfVariables);
X0 = [];
lb = -inf(numberOfVariables,1);
ub = inf(numberOfVariables,1);

defaultopt.ConvexCheck = 1;
defaultopt.MaxIter = 200;
defaultopt.TolFun = 1e-8;
defaultopt.TolX = 1e-12;
defaultopt.ConvexCheck = 1;

flags.computeLambda = computeFirstOrderOpt;
flags.detailedExitMsg = detailedExitMsg;
flags.verbosity = verbosity;
flags.caller = 'quadprog';

[X, fval] = ...
    ipqpdense(full(H), f, A, B, Aeq, Beq, lb, ub, X0, flags, ...
    options, defaultopt);
