% Estimate region of attraction by V-s-iteration and bisection.

import casos.toolboxes.sosopt.*
% system states
x = casos.PS('x',2,1);

% system dynamics
f = [x(2); -(1-x(1)^2)*x(2) - x(1)];


A = double(subs(nabla(f,x),x,zeros(2,1)));

P = lyap(A',eye(2));

Vinit = x'*P*x;

pcontour(Vinit,1,[-1 1 -1 1])
hold on
pcontour(Vinit,0.5,[-1 1 -1 1])

% Lyapunov function candidate
p = x'*x;
g = 3*x(2)^2 + x(1)^2 -1; 
% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,0));
s2 = casos.PS.sym('s2',monomials(x,2));

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

% options
opts = struct('sossol','mosek');

%% setup solver
% solver 1: gamma-step
sos1 = struct('x',[V; s1; s2; b],...
              'f',-b, ...
              'p',[]);

sos1.('g') = [s1; 
              s2; 
              V-l; 
              s2*(V-1)-nabla(V,x)*f-l; 
              s1*(p-b) + 1 - V];

% states + constraint are SOS cones
opts.Kx = struct('lin', 4);
opts.Kc = struct('sos', 5);
opts.verbose = 1;
    
Vlb  = casos.PS(basis(V),-inf);
Vub  = casos.PS(basis(V),+inf);
s1lb = casos.PS(basis(s1),-inf);
s1ub = casos.PS(basis(s1),+inf);
s2lb = casos.PS(basis(s2),-inf);
s2ub = casos.PS(basis(s2),+inf);
glb  = casos.PS(basis(b),-inf);
gub  = casos.PS(basis(b),+inf);

tic
    S1 = casos.nlsossol('S1','sequential',sos1,opts);
toc

tic
% solve
sol1 = S1('x0',[Vinit ; 1; x'*x ; 1], ...
          'lbx',[Vlb;s1lb;s2lb;glb], ...
          'ubx',[Vub;s1ub;s2ub;gub]);
toc

figure()
pcontour(g,0,[-1 1 -1 1],'k--')
hold on
pcontour(sol1.x(1),1,[-1 1 -1 1])

