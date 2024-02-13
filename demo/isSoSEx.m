% Polynomial optimization.

% indeterminate variable
x = casos.PS('x');
y = casos.PS('y');

% some polynomial
f = 2*x^4 + 2*x^3*y - x^2*y^2+ 5*y^4;

% define SOS problem:
%   min g s.t. (f + g) is SOS
sos = struct('x',[],'f',0,'g',f);

% constraint is scalar SOS cone
opts = struct('Kc',struct('s',1));

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);

% evaluate
sol = S();
switch (S.stats.UNIFIED_RETURN_STATUS)

    case 'SOLVER_RET_SUCCESS'
        disp('Polynomial is sos')
    otherwise
         disp('Could not find a SOS decomposition')
end

fmp = to_multipoly(f);

% Use SOSOPT to check if p is an SOS
opts = sosoptions();
[feas,z,Q,f] = issos(fmp,opts);