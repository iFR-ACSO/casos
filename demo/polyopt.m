% Polynomial optimization.
%% Sum-of-squares cone
% indeterminate variable
x = casos.Indeterminates('x');
% some polynomial
f = x^4 + 1*x;
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

fprintf('SOS: Minimum is %g.\n', full(sol.f))

%% Diagonally dominant sum-of-squares cone

% constraint is scalar DSOS cone
opts = struct('Kc',struct('dsos',1));

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);
% evaluate
sol = S();

fprintf('DSOS: Minimum is %g.\n', full(sol.f))

%% Scaled diagonally dominant sum-of-squares cone

% constraint is scalar SDSOS cone
opts = struct('Kc',struct('sdsos',1));

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);
% evaluate
sol = S();

fprintf('SDSOS: Minimum is %g.\n', full(sol.f))
