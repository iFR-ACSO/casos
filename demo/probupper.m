% Obtain upper bound on \int_K f(x) dx with
%   1. K a bounded semialgebraic set described by a single polynomial
%      inequality, k(x)>0
%   2. f(x) is a polynomial function which is nonnegative

% indeterminate variable
x = casos.Indeterminates('x', 2);

% nonnegative polynomial on K
f = x(1)*x(1) + x(2)^4 + x(2)^2 + x(1)*x(2);

% polynomial tha defines the bounded set K
k = x(1)^3-x(1)^4-x(2)^2; 

% K contained in B = [0 1] x [-0.4 0.4]
B = [0 1; -0.4 0.4];

% indicator-like function
v = casos.PS.sym('v',monomials(x,0:7),'gram');

% SOS multiplier
s = casos.PS.sym('q',monomials(x,0:10),'gram');

% obtain cost
cost = int(v*f,  x(1), B(1,1), B(1,2));
cost = int(cost, x(2), B(2,1), B(2,2));

% define SOS feasibility
sos = struct('x', [s; v], 'f', cost, 'g', -s*k+v-1);

% constraint is scalar SOS cone
opts = struct;
opts.Kx = struct('sos',2);
opts.Kc = struct('sos',1);

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);

% evaluate
sol = S();

% the result should be a value above 0.18162
fprintf('int_K f(x) dx = %d\n', full(sol.f)); 
