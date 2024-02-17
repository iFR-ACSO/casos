% Polynomial optimization.

% indeterminate variable
x = casos.PS('x');
% some polynomial
f = x^4 + 10*x;
% scalar decision variable
g = casos.PS.sym('g');

% define SOS problem:
%   min g s.t. (f + g) is SOS
sos = struct('x',g,'f',g,'g',f+g);
% constraint is scalar SOS cone
opts = struct('Kc',struct('s',1));
tic
% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);
% evaluate
sol = S();
toc
fprintf('Sedumi Minimum is %g.\n', double(sol.f))
tic
% solve by relaxation to SDP
S = casos.sossol('S','cdcs',sos,opts);
% evaluate
sol = S();

fprintf('CDCS Minimum is %g.\n', double(sol.f))
toc
