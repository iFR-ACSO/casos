% Estimate region of attraction by V-s-iteration and bisection.

import casos.toolboxes.sosopt.*
% system states
x = casos.PS('x',2,1);

% system dynamics
f = [x(2); -(1-x(1)^2)*x(2) - x(1)];


A = double(subs(nabla(f,x),x,zeros(2,1)));

P = lyap(A',eye(2));

Vinit = x'*P*x;



% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,0));
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,2));

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

% options
opts = struct('sossol','mosek');

g = 3*x(2)^2 + x(1)^2 -1; 

cost = dot(g - (V-b),g - (V-b)) - b^2;

%% setup solver
% solver 1: gamma-step
sos1 = struct('x',[V; s1; s2; s3;b],...
              'f',cost, ...
              'p',[]);

sos1.('g') = [s1; 
              s2; 
              s3;
              V-l; 
              s2*(V-b)-nabla(V,x)*f-l;
              s3*(V-b) - g];

% states + constraint are SOS cones
opts.Kx = struct('lin', 5);
opts.Kc = struct('sos', 6);
opts.verbose = 1;
    
Vlb  = casos.PS(basis(V),-inf);
Vub  = casos.PS(basis(V),+inf);
s1lb = casos.PS(basis(s1),-inf);
s1ub = casos.PS(basis(s1),+inf);
s2lb = casos.PS(basis(s2),-inf);
s2ub = casos.PS(basis(s2),+inf);
s3lb = casos.PS(basis(s3),-inf);
s3ub = casos.PS(basis(s3),+inf);
blb = casos.PS(basis(b),-inf);
bub = casos.PS(basis(b),+inf);

opts.Sequential_Algorithm = 'SQP';
tic
    S1 = casos.nlsossol('S1','sequential',sos1,opts);
toc

tic
% solve
sol1 = S1('x0',[Vinit^2 ; 1; x'*x; x'*x;1], ...
          'lbx',[Vlb;s1lb;s2lb;s3lb;blb], ...
          'ubx',[Vub;s1ub;s2ub;s3ub;bub]);
toc

figure()
pcontour(g,0,[-1 1 -1 1],'k--')
hold on
pcontour(sol1.x(1),double(sol1.x(end)),[-1 1 -1 1])

