% Obtain upperbound on vol(K) with K a semialgebraic set described by
% single polynomial inequality, k(x) > 0
% (use of stokes constraint to speed convergence)

% indeterminate variable
x = casos.Indeterminates('x', 2);

% polynomial tha defines the bounded set K
k = 1-x'*x;

% K contained in B = [0 1] x [-0.4 0.4]
B = [-1 1; -1 1];

% indicator-like function
v = casos.PS.sym('v',monomials(x,0:10),'gram');

% auxiliary polynomial for stokes constraint
w = casos.PS.sym('w',monomials(x,0:10));
divkw = sum(nabla(k*w,x));

% SOS multiplier
s = casos.PS.sym('q',monomials(x,0:10),'gram');

% obtain cost
cost = int(v,    x(1), B(1,1), B(1,2));
cost = int(cost, x(2), B(2,1), B(2,2));

% define SOS feasibility
sos = struct('x', [w; s; v], 'f', cost, 'g', -s*k+v+divkw-1);

% constraint is scalar SOS cone
opts = struct;
opts.Kx = struct('lin', 1, 'sos',2);
opts.Kc = struct('sos',1);

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);

% evaluate
sol = S();

% the result should be a value slightly above pi
fprintf('volume(K) = %d\n', full(sol.f)); 


