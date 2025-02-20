% Polynomial optimization.

% indeterminate variable
x = casos.Indeterminates('x');
% some polynomial
f = x^4 + 10*x;
% scalar decision variable
g = casos.PS.sym('g');

% define SOS problem:
%   min g s.t. (f + g) is SOS
sos = struct('x',g,'f',g,'g',f+g);
% constraint is scalar SOS cone
opts = struct('Kc',struct('sos',1));

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);
% evaluate
sol = S();

fprintf('Minimum is %g.\n', full(sol.f))
