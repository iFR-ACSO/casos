% Projection onto sum-of-squares cone

rng(0)

% indeterminate variable
x = casos.PS('x');
% some random polynomial
p = casos.PS(monomials(x,0:4),randn(5,1));

% Gram decision variable
s = casos.PS.sym('q',monomials(x,0:2),'gram');

% projection error
e = s - p;

% define Q-SOS problem:
%   min ||s-p||^2 s.t. s is SOS
sos = struct('x',s,'f',dot(e,e));
% states is scalar SOS cone
opts = struct('Kx',struct('sos',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
% evaluate
sol = S();

fprintf('Distance to SOS cone is %g.\n', double(sol.f))

