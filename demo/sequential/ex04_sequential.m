%% ------------------------------------------------------------------------
%
%
%   Short Description:  Calculate an inner-approximation of the 
%                       continous-time reachable set of the Van-der-Pol 
%                       Oscillator. Instead of a bisection we try to reduce
%                       the Euclidean distance to the scaled constraint set.
%                       Polynomial degrees and constraint similar to
%                       reference below.
%
%   Reference: Cunis, T., Kolmanovsky, I.- Viability, viscosity, and 
%              storage functions in model-predictive control with terminal 
%              constraints, Automatica, 2021, 
%              DOI: 10.1016/j.automatica.2021.109748
%           
%
%--------------------------------------------------------------------------

clear
close all
clc

import casos.toolboxes.sosopt.*


% system states
x = casos.PS('x',2,1);
u = casos.PS('u',1,1);
t  = casos.PS('t');
t0 = casos.PS.sym('t0');
t1 = casos.PS.sym('t1');

% system dynamics
f = [x(2);
        (1-x(1)^2)*x(2)-x(1)];

gx = [0;1];


% terminal region
P = [6.4314    0.4580
    0.4580    5.8227];

l = x'*P*x-1;              % l(x) \leq 0

% state constraint
g  = 3*x(2)^2 + x(1)^2 -1;  % g(x) \leq 0
g0 = 3*x(2)^4 + x(1)^4 -1;
% control constraint
umin = -1;
umax = 1;

% Time horizon4
T  = 1;

% time polynomial
hT = (t)*(T-t);

Vval = x'*P*x;

b  = casos.PS.sym('b');
V  = casos.PS.sym('v',monomials([x;t],0:5));
K  = casos.PS.sym('k',monomials([x;t],0:3));
s1 = casos.PS.sym('s1',monomials([x;t],2:4));
s2 = casos.PS.sym('s2',monomials([x;t],0:3));
s3 = casos.PS.sym('s3',monomials([x;t],0:3));
s4 = casos.PS.sym('s4',monomials([x;t],0:3));
s5 = casos.PS.sym('s5',monomials([x;t],0:3));
s6 = casos.PS.sym('s6',monomials([x;t],0:3));
s7 = casos.PS.sym('s7',monomials([x;t],0:3));
s8 = casos.PS.sym('s8',monomials([x;t],0:3));
s9 = casos.PS.sym('s9',monomials(x,0:3));
s10 = casos.PS.sym('s10',monomials(x,0:3));

% options
opts.sossol = 'mosek';
opts.max_iter = 100;
opts.verbose  = 1;

% zero sublevel set
b = 0;

sos = struct('x',[V; K; s1; s2; s3; s4; s5; s6;s7;s8;s9;s10],...
              'f',dot( g0-(subs(V,t,0)-b), g0-(subs(V,t,0)-b) )  , ...
              'p',[]);


% constraints
sos.('g') = [s1;
             s2; 
             s3;
             s4; 
             s5;
             s6;
             s7;
             s8;
             s9;
             s10;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*(f+gx*K);
             s3*(V-b) - s4*hT + K - umin;
             s5*(V-b) - s6*hT + umax - K;
             s7*(V-b) - s8*hT - g; 
             s9*l  - subs(V,t,0) + b;
             s10*(subs(V,t,T)-b)    - l 
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));

% build solver
solver_oneStepReach  = casos.nlsossol('S','filter-linesearch',sos,opts);

% initial guess
x0 = casos.PD([ Vval;  ...
                x'*x; ...
                x'*x; ...
                x'*x; ...
                x'*x;
                x'*x;
                x'*x;
                x'*x;
                x'*x;
                x'*x;
                x'*x;
                x'*x]);

% solve
sol = solver_oneStepReach('x0',x0); 

% augment solution with zero sublevel set (just for plotting)
sol.x = [sol.x;0];


%% plotting
import casos.toolboxes.sosopt.*

figure(1)
% reachability storage function at t = 0s
pcontour(subs(sol.x(1),t,0),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'b')
hold on 
% reachability storage function at T = 0s
pcontour(subs(sol.x(1),t,T),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'g--')

% constraint set
pcontour(g,0,[-1 1 -1 1],'k--')
% terminal set
pcontour(l,0,[-1 1 -1 1],'k')
legend('V(0,x)','V(T,x)','g(x)','l(x)')



