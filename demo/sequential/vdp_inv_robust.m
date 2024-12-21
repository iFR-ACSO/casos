
clear
close all
clc

% system states
x = casos.PS('x',2,1);
u = casos.PS('u',1,1);
w = casos.PS('w',1,1);

% to have enough margins
umin = -2;
umax =  2;

% dynamics
f = [
    x(2)
    (1-x(1)^2)*x(2) - x(1) + u
];

% trim point
x0    = [0 0]';
u0    = 0;

A = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));
   
[K0,P] = lqr(A, B, eye(2), 2.5);

f = [
    x(2)
    (1-x(1)^2)*x(2) - x(1) + u
] + [0;1]*w;


% %scale initial lyapunov
Vval = x'*P*x;

g = x(1)^2 + 3*x(2)^2 - 1;

V  = casos.PS.sym('v',monomials(x,2));
K  = casos.PS.sym('k',monomials(x,1),[length(u) 1]);

s1  = casos.PS.sym('s1', monomials([x;w],2));
s1r = casos.PS.sym('s1r',monomials([x;w],2));

s3  = casos.PS.sym('s3', monomials([x;w],0:2),[length(u) 1]);
s3r = casos.PS.sym('s3r',monomials([x;w],0:2),[length(u) 1]);

s5  = casos.PS.sym('s5', monomials([x;w],0:2),[length(u) 1]);
s5r = casos.PS.sym('s5r',monomials([x;w],0:2),[length(u) 1]);

s8  = casos.PS.sym('s8',monomials(x,0:2));

b  = casos.PS.sym('b');
% options
opts = struct('sossol','mosek');

% adjust optimality thresholds
% opts.conVioTol     = 6e-8;
% opts.optTol        = 1e-1;
% opts.error_on_fail = 0;
opts.verbose       = 1;
opts.max_iter      = 1000;

b = 0.5;
sos = struct('x',[V; K; s1; s3; s5;s8;s1r;s3r;s5r],...
              'f',dot(g-(V), ...
                      g-(V)) ,... 
              'p',[]);



wmax = 0.2;
% constraints
sos.('g') = [
             % s1;
             s3;
             s5;
             s8;
             s1r;s3r;s5r;
             s1*(V-b)  -  nabla(V,x)*subs(f,[x;u],[x;K]) + s1r*(w'*w-wmax^2);
             s3*(V-b)  + K - umin  + s3r*(w'*w-wmax^2)                      ; 
             s5*(V-b)  + umax - K  + s5r*(w'*w-wmax^2)                      ;
             s8*(V-b) - g                                 ;
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));

solver_oneStepReach  = casos.nlsossol('S','sequential',sos,opts);

x0 = casos.PD([ (Vval);(-K0*x);x'*x; x'*x; ...
                (x'*x); x'*x; x'*x; (x'*x);(x'*x)]);

sol = solver_oneStepReach('x0',x0);




%% check feasibility
Vsol    = sol.x(1)-b; 

%% plotting
import casos.toolboxes.sosopt.*


figure(1)
pcontour(Vsol,0,[-1 1 -1 1],'b')
hold on
pcontour(g,0,[-1 1 -1 1],'k--')
