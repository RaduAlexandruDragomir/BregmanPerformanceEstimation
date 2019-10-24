% In this script, we solve the semidefinite formulation of the Performance
% Estimation Problem (PEP) for the NoLips/Bregman Gradient algorithm.

clear all; clc;

%% Parameters of the problem
N       = 3;            % number of iterations 
L       = 1;            % h-smoothness constant
lambda  = 0.8/L;    % step size
R       = 1;            % initial radius

% solver parameters
verbose     = 1;
tolerance   = 1e-8;
tolerance2   = 1e-4;

%% Setting up the problem variables

% P = [ x0 ... xN | g0 ... gN | s0]
% G = P^T P
% F = [             f0 ... fN ]
% H = [             h0 ... hN ]

dimG  = 2*N + 3;
dimF  = N + 1;
dimH  = dimF;
nbPts = N + 2; % x*, x0, ...,xN

% encoding vectors
% for example, x(i,:)' * G * x(j,:) encodes <xi, xj>
x = zeros(nbPts, dimG);
g = zeros(nbPts, dimG);
s = zeros(nbPts, dimG);
f = zeros(nbPts, dimF);
h = zeros(nbPts, dimH);

x(2:nbPts,1:nbPts-1) = eye(nbPts-1);
g(2:nbPts,N+2:2*N+2) = eye(N+1);
s(2,2*N+3)           = 1;


f(2:nbPts, 1:nbPts-1)      = eye(nbPts-1);
h                          = f;

% encoding the NoLips algorithm
for i = 1:N
    s(2+i,:) = s(1+i,:) - lambda * g(1+i,:);
end    

% for convenience, define encoding vectors corresponding to the optimum x*
xs = x(1,:); gs = g(1,:); ss = s(1,:); fs = f(1,:); hs = h(1,:);

% and encoding vectors corresponding to the points x0...xN
xk = x(2:end,:); gk = g(2:end,:); sk = s(2:end,:); fk = f(2:end,:); hk = h(2:end,:);

%% setting up the SDP

G = sdpvar(dimG);       % G is dimG x dimG symmetric
constraint = ( G >= 0); % G is PSD
F = sdpvar(dimF,1);     % F is  dimF x1
H = sdpvar(dimF,1);     % H is  dimF x1

% initial radius constraint
constraint = constraint + ( L*(hs-hk(1,:))*H - L*sk(1,:)*G*(xs-xk(1,:))'<= R);

% convexity constraints
for i = 1:nbPts         % must be used for xs,x0,...,xN 
    for j = 1:nbPts     % must be used for xs,x0,...,xN
        if i ~= j
            % convexity of f
            constraint = constraint + ((f(j,:)-f(i,:))*F + g(j,:)*G*(x(i,:)-x(j,:))'  <= 0);
            
            % convexity of Lh - f
            constraint = constraint + ( -(L*(h(i,:)-h(j,:))*H + (f(j,:)-f(i,:))*F - (L*s(j,:)-g(j,:))*G*(x(i,:)-x(j,:))') <= 0);
        end
    end
end

% objective value
obj = (fk(end,:)-fs)*F;

% solving the PEP
solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(constraint,-obj,solver_opt);

fprintf("\nProblem info: %s\n\n", solverDetails.info)
fprintf("PEP value:          %d\n", double(obj))
fprintf("theoretical value:  %d\n", R / lambda / N)




