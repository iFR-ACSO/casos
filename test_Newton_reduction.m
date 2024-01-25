%% test Newton polytope monomial basis reduction reduction
% Go to grambasis.m to uncomment the code to enable the newton polytope 

x = casos.PS('x', 3, 1);
f = 12+x(2)^2-2*x(1)^3*x(2)+2*x(2)*x(3)^2+x(1)^6-2*x(1)^3*x(3)^2+x(3)^4+x(1)^2*x(2)^2;
[Z,K,z] = grambasis(f);


%% Verify if f is SOS
sos = struct('x', [] ,'g',f);
% states + constraint are SOS cones
opts.Kc.s = 1;

% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);
tic
sol = S();
ftime = ftime + toc;

fprintf('%s: %ds \n', S.stats.UNIFIED_RETURN_STATUS, ftime)

