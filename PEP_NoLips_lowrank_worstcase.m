% Script for examining the worst-case function for the NoLips/Bregman gradient
% algorithm. We use the PEP formulation along with a rank minimization heuristic.

clear all; clc;

%% Parameters of the problem
% Note that the example will be displayed in a plot only 
% for N = 1,2 (because it will be in dimension 2 and 3 respectively)
N       = 1;            % number of iterations 
L       = 1;            % h-smoothness constant
lambda  = 1/L;    % step size
R       = 1;            % initial radius

% solver parameters
verbose     = 1;
tolerance   = 1e-8;

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

% the algorithm

for i = 1:N
    s(2+i,:) = s(1+i,:) - gamma(i) * g(1+i,:);
end


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
constraint = constraint + ( (hs-hk(1,:))*H - sk(1,:)*G*(xs-xk(1,:))'<= R);

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

% we now search for a "good" solution. First add the constraint
% that we have a maximizer of the problem. We can do this because
% we know that the problem value is r / lambda / N
constraint = constraint + ((fk(end,:)-fs)*F >= R / lambda /N);

% then add some orthogonality constraints on the gradients of f
% to get a "worst function in the world" type behavior
for i = 2:N+1
    for j = i+1:N+2
       constraint = constraint + (g(i,:)*G*(g(j,:))' == 0);
    end
    
   
end

% we also can impose this normalization constraint in order to get
% "  coordinates " (falcutative)
constraint = constraint + (x(2,:)*G*x(2,:)' == N+1);
constraint = constraint + (x(3,:)*G*x(3,:)' == N);

% we search for the minimal trace solution, hopefully leading to a low rank
obj = trace(G) ;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(constraint, obj,solver_opt);

%% we perform QR decomposition to recover P from G = P' * P 
[V,D]           =eig(double(G));
tol_eigvalues   =1e-4; %Throw away eigenvalues smaller that tol
eigenV          =diag(D); 
eigenV(eigenV < tol_eigvalues)   =0;
new_D=diag(eigenV); [~,P]=qr(sqrt(new_D)*V.');

P=P(1:sum(eigenV>0),:);

dim = size(P,1);

% for a more convenient representation, we rotate P
% on the basis given by grad d(x1)
gr_d = (P * (L*s' - g'));
grad_dx1 = gr_d(:,2);
u1 = grad_dx1 / norm(grad_dx1);

if dim == 2
    
    rotation = [u1(2) -u1(1);
                u1(1) u1(2)];     
            
    P = rotation * P;
end
        
% we can now recover the problem data thanks to P and the encoding vectors
X  = P * x'; % x_i
Gf = P * g'; % grad f(x_i)
Gh = P * s'; % grad h(x_i)

Gd = L * Gh - Gf; % grad d(x_i)

% if the dimension is 1 or 2, we can plot the functions
if dim == 1 || dim == 2
    plot_discrete_functions(X,Gf,Gh,F,H)
end

