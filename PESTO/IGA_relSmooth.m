clear all; clc;

% In this example, we try to certify the improved interior gradient
% algorithm (IGA) from [1]use the improved interior gradient algorithm (IGA) 
% in the setting where f is convex and smooth relative to a kernel h
%
% In IGA.m, we certify this method when f is assumed to be L-smooth and
% h strongly convex. Here, we try the more general h-smooth assumption.
% We show that any value of t_k that is great than 1 causes the PEP to
% be unbounded, hence indicating that IGA fails converge.
% 
% The method originates from:
% [1] Alfred Auslender, and Marc Teboulle. "Interior gradient and proximal
%     methods for convex and conic optimization."
%     SIAM Journal on Optimization (2006).
%

% (0) Initialize an empty PEP
P = pep();


% (1) Set up the objective function
L = 1; % relative smoothness constant
d  = P.DeclareFunction('Convex'); % d = Lh - f1 is convex
f1 = P.DeclareFunction('Convex');
h  = (d + f1)/L;

F  = f1;

% (2) Set up the starting point and initial condition
x0        = P.StartingPoint();     % x0 is some starting point
[xs,fs]   = F.OptimalPoint('opt'); % xs is an optimal point, and fs=F(xs)
[sxs,hxs] = h.oracle(xs, 'opt');                  
[sx0,hx0] = h.oracle(x0, 'x0');                  
[gx0,fx0] = f1.oracle(x0, 'x0');

% (3) Algorithm
lambda  = 1/L;          % stepsize
N       = 3;           % number of iterations

% setting up the step size growth: any value that is higher than 1
% will cause the PEP to be unbounded, indicating that the algorithm
% diverges
tk = @(k) 1.1;

% store iterates in cells:
x  = cell(N+1,1);  x{1} = x0; 
z  = cell(N+1,1);  z{1} = x0;
y  = cell(N,1); 
g  = cell(N,1); g{1} = gx0;    % gradients of f at the y's
f  = cell(N,1); f{1} = fx0;    % function values of f at the y's
sz = cell(N+1,1); sz{1} = sx0; % subgradients of h at the z's 

for i=1:N
    
    y{i}   = (1- 1/tk(i) ) * x{i} + 1 / tk(i) * z{i};

    if i>1 % when i = 1, y{i} = x{i} so evaluation is useless
        [g{i},f{i}] = f1.oracle(y{i});
    end
    name    = sprintf('z%d',i);
    z{i+1}  = mirror(g{i}, sz{i}, h, lambda * tk(i), name); 
    x{i+1} = (1- 1/tk(i) ) * x{i} + 1 / tk(i) * z{i+1};
    
    [sz{i+1}, ~ ] = h.oracle(name); % this is for the next mirror step
end

Dh0 = hxs - hx0 - sx0 * (xs-x0);
P.InitialCondition(Dh0 + fx0-fs<=1); % Add an initial condition 

% (4) Set up the performance measure
fN = F.value(x{N+1});         % fN=F(xN)
P.PerformanceMetric(fN-fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(fN-fs)]  % worst-case objective function accuracy