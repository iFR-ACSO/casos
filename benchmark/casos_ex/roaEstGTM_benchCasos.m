%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% CaSoS. The quasi-convex solver is used to perform the bisections.
%
%--------------------------------------------------------------------------

function [gval,bval, solverTime_total,buildTime,callTime] = roaEstGTM_benchCasos()

%% problem definition
% system states 
x = casos.Indeterminates('x',4);

% Polynomial Dynamics
f = GTM_dynamics(x(1),x(2),x(3),x(4));

% shape function
p = x'*x*1e2;

% initial Lyapunov function candidate
P = [395.382671347059	-23.0032507976836	3.16965275615691	29.2992065909380
    -23.0032507976836	90.4764915483638	-16.1191428789579	-132.594376986429
    3.16965275615691	-16.1191428789579	3.44002648214490	24.2292666551058
    29.2992065909380	-132.594376986429	24.2292666551058	202.114797577027];

Vval = x'*P*x;

% enforce positivity
l = 1e-6*(x'*x);

% start to measure build time of all parameterized solver
buildTime_start = tic;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,1:2),'gram');
s2 = casos.PS.sym('s2',monomials(x,0:2),'gram');

% level of stability
g = casos.PS.sym('g');
b = casos.PS.sym('b');

% options
opts               = struct('sossol','mosek');
opts.error_on_fail = 0;
opts.verbose       = 0;
opts.conf_interval = [-1000 0];

%% setup solver
% solver 1: gamma-step

sos1 = struct('x',s1, ... % dec.var
              'f',-g, ... % cost function for bisection
              'p',V);     % parameter

% constraint
sos1.('g') = s1*(V-g)-nabla(V,x)*f-l;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

% build first solver
S1 = casos.qcsossol('S1','bisection',sos1,opts);

% solver 2: beta-step
sos2 = struct('x',s2, ... % dec.var
              'f',-b, ... % cost function for bisection
              'p',[V;g]); % parameter

% constraint
sos2.('g') = s2*(p-b)+g-V;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

% build second solver
S2 = casos.qcsossol('S2','bisection',sos2,opts);

% solver 3: V-step
sos3 = struct('x',V, ...        % dec.var
              'p',[b,g,s1,s2]); % parameter

% constraints
sos3.('g') = [V-l; 
              s2*(p-b)+g-V; 
              s1*(V-g)-nabla(V,x)*f-l];

opts    = struct;
opts.Kx = struct('sos', 0, 'lin', 1); 
opts.Kc = struct('sos', 3);

% build third solver
S3 = casos.sossol('S','mosek',sos3,opts);

% measure net build time
tmpbuildTime = toc(buildTime_start);

% initialize arrays
solvetime_all1 = zeros(100,1);
solvetime_all2 = zeros(100,1);
solvetime_all3 = zeros(100,1);

callTime1 = zeros(100,1);
callTime2 = zeros(100,1);
callTime3 = zeros(100,1);

% needed for convergence check
bval_old = [];

%% V-s-iteration
for iter = 1:100

    %% gamma step
    startSolve1 = tic;
    
    % call solver for gamma-step
    sol1 = S1('p',Vval);    

    % get all low level solver times from all subiterations of the biscetion
    solvetime_all1(iter) = sum(cellfun(@(x) x.solvetime_matlab, S1.stats.iter));
    
    % measure the CaSoS solver time i.e. call time --> additional "build" time
    callTime1(iter)     = toc(startSolve1) - solvetime_all1(iter);

    % extract solution
    gval = -sol1.f;
    s1val = sol1.x;

    %% beta step
    startSolve2 = tic;

    % call solver beta-step
    sol2 = S2('p',[Vval;gval]); 

    % get all low level solver times from all subiterations of the biscetion
    solvetime_all2(iter) = sum(cellfun(@(x) x.solvetime_matlab, S2.stats.iter));

    % measure the CaSoS solver time i.e. call time --> additional "build" time
    callTime2(iter)     = toc(startSolve2) - solvetime_all2(iter);
    
    % extract solution
    bval = -sol2.f;
    s2val = sol2.x;

    %% V-step
    startSolve3 = tic;

    % call solver V-step
    sol3 = S3('p',[bval,gval,s1val,s2val]);

    % get low level solver time
    solvetime_all3(iter) = S3.stats.solvetime_matlab;
    
    % measure the CaSoS solver time i.e. call time --> additional "build" time
    callTime3(iter) = toc(startSolve3) - solvetime_all3(iter);

    % extract solution
    Vval = sol3.x;

    % show progress 
    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(gval));
    
    %% check convergence via shape function size
    if ~isempty(bval_old)
        if abs(full(bval-bval_old)) <= 1e-3
            break
        else
            bval_old = bval;
        end
    else
        bval_old = bval;
    end

end % end for-loop

%% prepare output

callTime  = sum(callTime1) + sum(callTime2) + sum(callTime3);
buildTime = tmpbuildTime + callTime;

% total solver time over all iterations
solverTime_total = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);

% save the complete workspace, so people do not have to re-run execpt they
% want to
% save('Casos_GTM_ROA_bench.mat')

end % end of function

