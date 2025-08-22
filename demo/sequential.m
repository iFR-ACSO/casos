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

close all
clear
clc

import casos.toolboxes.sosopt.*

%% system states
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

% minimize the quadratic distance to a given sublevel set
p = x'*x;

%% setup solver

% options
opts = struct('sossol','mosek');
opts.verbose = 1;
opts.hessian_approx = 'BFGS';
opts.Hessian_init   = 'Identity';
opt.scale_BFGS0     = 1e-3;

sos = struct('x',[V; s2;s1;b],...
              'f',-b, ...
              'p',[]);

% constraints
sos.('g') = [s2; 
             s1;
             V-l; 
             s2*(V-1)-nabla(V,x)*f-l;
             s1*(p-b) + 1 - V];



% states + constraint are linear/SOS cones
opts.Kx = struct('lin', 4);
opts.Kc = struct('sos', 5);

% build sequential solver
solver  = casos.nlsossol('S','sequential',sos,opts);


%% solve
sol = solver('x0',[ Vinit;  (x'*x);  1; 1]); 

%% plot sublevel set
figure()
Vval = sol.x(1);
beta_val = full(casadi.DM(full(sol.x(end))));
pcontour(sol.x(1),1,[-4 4 -4 4])
hold on
pcontour(p,beta_val,[-4 4 -4 4],'r')


