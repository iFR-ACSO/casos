%% ------------------------------------------------------------------------
%
%
%   Short Description:  Demo for sequential SOS using a line search filter. 
%                       Calculate an inner-estimate of the
%                       region-of-attraction for the closed-loop
%                       Van-der-Pol Oscilator.
%
%   Reference: -
%           
%
%--------------------------------------------------------------------------

% system states
x = casos.Indeterminates('x',2,1);

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% initial Lyapunov function
Vinit = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials(x,2));
s1 = casos.PS.sym('s1',monomials(x,0));
% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

% quadratic shape function
p = x'*x;

%% Setup solver
% options
opts = struct('sossol','mosek');
opts.verbose = 1;
opts.hessian_approx = 'BFGS';
opts.hessian_init   = 'Identity';
opts.scale_BFGS0     = 1e-3;

sos = struct('x',[V;s2;s1;b],'f',-b,'p',[]);

% constraints
sos.('g') = [
             s2;
             s1;
             V-l; 
             s2*(V-1)-nabla(V,x)*f-l;
             s1*(p-b) + 1 - V
];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', 4);
opts.Kc = struct('sos', 5);

% build sequential solver
S = casos.nlsossol('S','sequential',sos,opts);


%% Solve nonlinear problem
sol = S('x0',[Vinit; (x'*x); 1; 1]); 

%% Plot results
% Lyapunov function
Vfun = to_function(sol.x(1));
% shape function
pfun = to_function(p);
% stable level set
bsol = full(sol.x(end));

figure
fcontour(@(x,y) full(Vfun(x,y)), [-2 2], 'b', "LevelList", [1 1])
hold on
fcontour(@(x,y) full(pfun(x,y)), [-2 2], 'r--', "LevelList", [bsol bsol])
hold off
grid on
