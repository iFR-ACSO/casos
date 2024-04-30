% Estimate region of attraction by V-s-iteration and bisection.

import casos.toolboxes.sosopt.*
% system states
x = casos.PS('x',2,1);

% system dynamics
f = [x(2); -(1-x(1)^2)*x(2) - x(1)];


A = double(subs(nabla(f,x),x,zeros(2,1)));

P = lyap(A',0.1*eye(2));

Vinit = x'*P*x;


% Lyapunov function candidate
V = casos.PS.sym('v',basis(cleanpoly(sol1.x(1),1e-5)));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,2));
s2 = casos.PS.sym('s2',basis(cleanpoly(sol1.x(1),1e-5)));
s3 = casos.PS.sym('s3',monomials(x,2));

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

% options
opts = struct('sossol','mosek');

g = 3*x(2)^2 + x(1)^2 -1; 

cost = dot(g - (V-1),g - (V-1)) ;

%% setup solver
% solver 1: gamma-step
sos1 = struct('x',[V; s2; s3],...
              'f',cost, ...
              'p',[]);

sos1.('g') = [s2; 
              s3;
              V-l; 
              s2*(V-1)-nabla(V,x)*f-l;
              s3*(V-1) - g];

% states + constraint are SOS cones
opts.Kx = struct('lin', 3);
opts.Kc = struct('sos', 5);
opts.verbose = 1;
    
Vlb  = casos.PS(basis(V),-inf);
Vub  = casos.PS(basis(V),+inf);
s1lb = casos.PS(basis(s1),-inf);
s1ub = casos.PS(basis(s1),+inf);
s2lb = casos.PS(basis(s2),-inf);
s2ub = casos.PS(basis(s2),+inf);
s3lb = casos.PS(basis(s3),-inf);
s3ub = casos.PS(basis(s3),+inf);
blb  = casos.PS(basis(b),-inf);
bub  = casos.PS(basis(b),+inf);

opts.Sequential_Algorithm = 'SQP';
tic
    S1 = casos.nlsossol('S1','sequential',sos1,opts);
toc

tic
% solve
sol1 = S1('x0',[Vinit ; x'*x; x'*x], ...
          'lbx',[Vlb;s2lb;s3lb], ...
          'ubx',[Vub;s2ub;s3ub]);
toc

figure()
pcontour(g,0,[-1 1 -1 1],'k--')
hold on
pcontour(sol1.x(1),1,[-1 1 -1 1])

