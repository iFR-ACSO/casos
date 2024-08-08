% Projection onto sum-of-squares cone

rng(0)

% indeterminate variable
x = casos.Indeterminates('x');

% some random polynomial
s1 = casos.PD(monomials(x,0:4),randn(5,1));
V = casos.PD(monomials(x,0:4),randn(5,1));
p = casos.PD(monomials(x,0:4),randn(5,1));

% constraint
c = s1*V-p;

% Gram decision variable
s = casos.PS.sym('q',sparsity(c));

% projection error
e = s - p;

% define Q-SOS problem:
%   min ||s-p||^2 s.t. s is SOS

sos = struct('x',s,'f',dot(e,e),'g',s);

% states is scalar SOS cone
opts = struct('Kx',struct('lin',1),'Kc',struct('sos',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);

tic
% evaluate
sol = S();
toc

fprintf('Distance to SOS cone is %g.\n', full(sol.f))

%% now only gram
opts = [];
% Gram decision variable
s = casos.PS.sym('q',grambasis(sparsity(c)));

% projection error
e = s - c;

% define Q-SOS problem:
%   min ||s-p||^2 s.t. s is SOS
sos = struct('x',s,'f',dot(e,e));

% states is scalar SOS cone
opts = struct('Kx',struct('sos',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
tic
% evaluate
sol = S();
toc
fprintf('Distance to SOS cone is %g.\n', full(sol.f))

