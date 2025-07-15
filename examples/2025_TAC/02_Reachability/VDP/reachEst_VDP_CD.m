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
umax =  1;

% Time horizon
T  = 1;

% time polynomial
hT = (t)*(T-t);


b  = casos.PS.sym('b');

V  = casos.PS.sym('v',monomials([x;t],0:2));
K  = casos.PS.sym('k',monomials([x;t],0:2));
s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],0:2));
s3 = casos.PS.sym('s3',monomials([x;t],0:2));
s4 = casos.PS.sym('s4',monomials([x;t],0:2));
s5 = casos.PS.sym('s5',monomials([x;t],0:2));
s6 = casos.PS.sym('s6',monomials([x;t],0:2));
s7 = casos.PS.sym('s7',monomials([x;t],0:2));
s8 = casos.PS.sym('s8',monomials([x;t],0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));
s10 = casos.PS.sym('s10',monomials(x,0:2));


% opts = struct('sossol','mosek');
% opts.verbose = 1;
% opts.error_on_fail = 0;
% opts.conf_interval = [0 1];
% zero sublevel set

sos = struct('x',[K; s1; s2; s3; s4; s5; s6;s7;s8;s9],...
              'p',[V;b]);


% constraints
sos.('g') = [
             s1;
             s2; 
             s3;
             s4; 
             s5;
             s6;
             s7;
             s8;
             s9;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*(f+gx*K);
             s3*(V-b) - s4*hT + K - umin;
             s5*(V-b) - s6*hT + umax - K;
             s7*(V-b) - s8*hT - g; 
             s9*(subs(V,t,T)-b)    - l;
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));

% build solver
solver_Reach1  = casos.sossol('S','mosek',sos,opts);

opts = [];
sos = struct('x',[V;s10;b],...
              'f', dot( g0-(subs(V,t,0)-b), g0-(subs(V,t,0)-b) ) , ...
              'p',[K; s1; s2; s3; s4; s5; s6;s7;s8;s9]);


% constraints
sos.('g') = [
             s10
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*(f+gx*K);
             s3*(V-b) - s4*hT + K - umin;
             s5*(V-b) - s6*hT + umax - K;
             s7*(V-b) - s8*hT - g; 
             s9*(subs(V,t,T)-b)    - l;
             s10*l  - subs(V,t,0) + b;
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));

solver_Reach2  = casos.sossol('S','mosek',sos,opts);



Vval = l;
bval = 0.01;
fval_old = [];
for iter = 1:100

% solve
sol1 = solver_Reach1('p',[Vval;bval]); 


sol2 = solver_Reach2('p',sol1.x); 

Vval = sol2.x(1);
bval = sol2.x(3);


% check convergence of cost function
 if ~isempty(fval_old)
    if abs(full(sol2.f-fval_old)) <= 1e-3
        break
    else
        fval_old = sol2.f;
    end
else
    fval_old = sol2.f;
 end

 % show progress 
fprintf('Iteration %d: f = %g, b = %g.\n',iter,full(sol2.f),full(bval));

end

%% plotting
import casos.toolboxes.sosopt.*

figure(1)
% reachability storage function at t = 0s
pcontour(subs(sol2.x(1),t,0),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'b')
hold on 
% reachability storage function at T = 0s
pcontour(subs(sol.x(1),t,T),full(casadi.DM(full(sol.x(end)))),[-1 1 -1 1],'g--')

% constraint set
pcontour(g,0,[-1 1 -1 1],'k--')
% terminal set
pcontour(l,0,[-1 1 -1 1],'k')
legend('V(0,x)','V(T,x)','g(x)','l(x)')



