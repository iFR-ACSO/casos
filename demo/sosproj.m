% Projection onto sum-of-squares cone

rng(0)

% indeterminate variable
x = casos.Indeterminates('x');
% some random polynomial
mono = monomials(x,2);
p = casos.PD(mono,randn(mono.nnz,1));

% Gram decision variable
s = casos.PS.sym('q',grambasis(p));

% projection error
e = s - p;

% define Q-SOS problem:
%   min ||s-p||^2 s.t. s is SOS
sos = struct('x',s, ...
             'f',dot(e,e), ...
             'g',s);
% states is scalar SOS cone
opts = struct('Kx',struct('lin',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
% evaluate
sol = S();

fprintf('Distance to SOS cone is %g.\n', full(sol.f))
