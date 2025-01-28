%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% 
%
%
%

% Inner-approximation reachable set Van-der-Pol Oscillator

clear
close all
clc

import casos.toolboxes.sosopt.*

profile on
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
% opts.hessian_approx = 'BFGS';
% adjust optimality thresholds
% opts.conVioTol = 1;
% opts.optTol    = 1e-1;
% opts.error_on_fail = 0;
opts.verbose = 1;
% opts.conViolCheck = 'pseudo';
opts.max_iter = 100;
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
% profile on
% build sequential solver
buildTime_in = tic;
    solver_oneStepReach  = casos.nlsossol('S','filter-linesearch',sos,opts);
buildtime = toc(buildTime_in);
% profile viewer
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

sol = solver_oneStepReach('x0',x0); 

sol.x = [sol.x;0];

% sol = solver_oneStepReach('x0',sol.x); 

disp(['Solver buildtime: ' num2str(buildtime), ' s'])


%%
% isSOS(sol.g(1))
% isSOS(sol.g(2))
% isSOS(sol.g(3))
% isSOS(sol.g(4))
% isSOS(sol.g(5))
% isSOS(sol.g(6))
% isSOS(sol.g(7))
% isSOS(sol.g(8))
% % isSOS(sol.g(9))
% isSOS(sol.g(10))
% isSOS(sol.g(11))
% isSOS(sol.g(12))
% isSOS(sol.g(13))
% isSOS(sol.g(14))
% isSOS(sol.g(15))
profile viewer

%% plotting
import casos.toolboxes.sosopt.*

figure(1)
pcontour(subs(sol.x(1),t,0),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'b')
hold on 
pcontour(g,0,[-1 1 -1 1],'k--')
pcontour(l,0,[-1 1 -1 1],'k')
pcontour(subs(sol.x(1),t,T),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'g--')

% for dt = 0.1:0.1:1.9
%     pcontour(subs(sol.x(1),t,dt),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'m')
% end


