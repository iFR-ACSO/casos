% verify if polynomial is sos
x = casos.PS('x',3,1);
p = 12+x(2)^2-2*x(1)^3*x(2)+2*x(2)*x(3)^2+x(1)^6-2*x(1)^3*x(3)^2 ...
    +x(3)^4+x(1)^2*x(2)^2;

% define sos feasibility
sos = struct('x',[],'g',p);
opts.Kc.s = 1;
opts.Kx.s = 0;

% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);

% evaluate parametrized SOS problem
tic
sol = S();
toc
disp(S.stats.UNIFIED_RETURN_STATUS)