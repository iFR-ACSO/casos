%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% CaSoS. The quasi-convex solver is used to perform the bisections.
%
%--------------------------------------------------------------------------

function [solverTimes,buildTimes,gval_array] = roaEstNlinkPend_benchVsIterCasos(Nmax,deg)


buildTimes  = zeros(Nmax-1,1);
solverTimes = zeros(Nmax-1,1);
gval_array  = zeros(Nmax-1,1);

for n = 2:Nmax
    % system states
    x = casos.PS('x',2*n,1);
    
    % system dynamics
    f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);
    load(['data_n' num2str(n)])

    % Lyapunov function candidate
    Vval = x'*S*x;
    p = Vval*10;
    
  % enforce positivity
l = 1e-6*(x'*x);


% start to measure build time of all parameterized solver
buildTime_start = tic;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,1:2),'gram');
s2 = casos.PS.sym('s2',monomials(x,0:2),'gram');


% level of stability
b = casos.PS.sym('b');



%% setup solver
% solver 1: gamma-step

sos1 = struct('x',s1, ... % dec.var
              'p',V);     % parameter

% constraint
sos1.('g') = s1*(V-1)-nabla(V,x)*f-l;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);
% 
% % build first solver
S1 = casos.sossol('S1','mosek',sos1,opts);


opts               = struct('sossol','mosek');
opts.error_on_fail = 0;
opts.conf_interval = [-1000 0];
% solver 2: beta-step
sos2 = struct('x',s2, ... % dec.var
              'f',-b,...
              'p',V); % parameter

% constraint
sos2.('g') = s2*(p-b)+1-V;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

% build second solver
S2 = casos.qcsossol('S2','bisection',sos2,opts);

% solver 3: V-step
sos3 = struct('x',V, ...        % dec.var
              'p',[b,s1,s2]); % parameter

% constraints
sos3.('g') = [V-l; 
              s2*(p-b)+1-V; 
              s1*(V-1)-nabla(V,x)*f-l];

opts     = struct;
opts.Kx = struct('sos', 0, 'lin', 1); 
opts.Kc = struct('sos', 3);

% build third solver
S3 = casos.sossol('S','mosek',sos3,opts);

buildTimes(n-1) = toc(buildTime_start);

% initialize arrays
solvetime_all1 = zeros(100,1);
solvetime_all2 = zeros(100,1);
solvetime_all3 = zeros(100,1);

bval_old = [];

%% V-s-iteration
disp(['Start V-s-iteration for ' num2str(n) '-link pendulum'])
for iter = 1:20

    % gamma step
    sol1 = S1('p',Vval);    

    % get all solver times from all subiterations of the biscetion
    solvetime_all1(iter) = S1.stats.solvetime_matlab;

    % extract solution
    s1val = sol1.x;

    % beta step
    sol2 = S2('p',Vval); 

    % get all solver times from all subiterations of the biscetion
    solvetime_all2(iter) = sum(cellfun(@(x) x.solvetime_matlab, S2.stats.iter));

    % extract solution
    bval = -sol2.f;
    s2val = sol2.x;

    % V-step
    sol3 = S3('p',[bval,s1val,s2val]);
    
    % extract solution
    Vval = sol3.x;

    % get solver time
    solvetime_all3(iter) = S3.stats.solvetime_matlab;
    
    % show progress 
    fprintf('Iteration %d: b = %g.\n',iter,full(bval));
    
    % check convergence
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

% total solver time over all iterations
solverTimes(n-1) = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);
end

end % end of function