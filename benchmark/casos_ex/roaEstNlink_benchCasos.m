%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the N-link pendulum dynamics
% in CaSoS.
%
%--------------------------------------------------------------------------

function [bval_array, solverTimes_total,buildTimes] = roaEstNlink_benchCasos(Nmax,deg)

buildTimes          = zeros(Nmax-1,1);
solverTimes_total   = zeros(Nmax-1,1);
bval_array          = zeros(Nmax-1,1);

% iterate of the Number of links
for n = 2:Nmax

disp(['Compute maximum ROA for the ' num2str(n) '-link pendulum'])

x = casos.PS('x',2*n,1);

% get dynamics
f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);

% load P matrix
load(['data_n' num2str(n)])

% Lyapunov function candidate
Vval = x'*P*x;
p = x'*x*100;
    
% enforce positivity
l = 1e-6*(x'*x);

% start to measure build time of all parameterized solver
buildTimes_start = tic;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,1:2),'gram');
s2 = casos.PS.sym('s2',monomials(x,0:1),'gram');

% level of stability
b = casos.PS.sym('b');

% fixed gamma
gval = 1;

%% setup solver

% solver 1: s1-step
sos1 = struct('x',s1,'p',V);     % parameter

% constraint
sos1.('g') = s1*(V-gval)-nabla(V,x)*f-l;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

% build first solver
S1 = casos.sossol('S1','mosek',sos1,opts);

% options
opts               = struct('sossol','mosek');
opts.error_on_fail = 0;
opts.verbose       = 0;

% solver 2: beta-step
sos2 = struct('x',s2, ... % dec.var
              'f',-b, ... % cost function for bisection
              'p',V); % parameter

% constraint
sos2.('g') = s2*(p-b)+gval-V;

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
              s2*(p-b)+gval-V; 
              s1*(V-gval)-nabla(V,x)*f-l];

opts     = struct;
opts.Kx = struct('sos', 0, 'lin', 1); 
opts.Kc = struct('sos', 3);

% build third solver
S3 = casos.sossol('S','mosek',sos3,opts);

tempBuildTime = toc(buildTimes_start);

% initialize arrays
solvetime_all1 = zeros(100,1);
solvetime_all2 = zeros(100,1);
solvetime_all3 = zeros(100,1);

callTime1 = zeros(100,1);
callTime2 = zeros(100,1);
callTime3 = zeros(100,1);

bval_old = [];

%% V-s-iteration
for iter = 1:20

    %% s1-step
    s1start = tic;
     
    % call solver
    sol1 = S1('p',Vval);    

    solvetime_all1(iter) = S1.stats.solvetime_matlab;
    callTime1(iter)      = toc(s1start)-solvetime_all1(iter);

    % extract solution
    s1val = sol1.x;

    %% beta step
    s2start = tic;

    % call solver
    sol2 = S2('p',Vval); 

    % get all solver times from all subiterations of the biscetion
    solvetime_all2(iter) = sum(cellfun(@(x) x.solvetime_matlab, S2.stats.iter));
    callTime2(iter)      = toc(s2start)-solvetime_all2(iter);

    % extract solution
    bval = -sol2.f;
    s2val = sol2.x;

    %% V-step
    s3start = tic;

    sol3 = S3('p',[bval,s1val,s2val]);

    solvetime_all3(iter) = S3.stats.solvetime_matlab;
    callTime3(iter) = toc(s3start)-solvetime_all3(iter);

    % extract solution
    Vval = sol3.x;

    % show progress 
    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(1));
  

end % end for-loop
 
% import casos.toolboxes.to_multipoly
% % we use multipoly function to compute the volume also for CaSoS
% V_c_mult =  to_multipoly(Vval-1);
% p_mult =  to_multipoly(p-bval);
% figure(3)
% pvar x_1 x_2 x_3 x_4 t % for plotting
% pcontour(subs(V_c_mult,[x_3;x_4],[0;0]),0,[-1 1 -1 1],'r')
% hold on
% pcontour(subs(p_mult,[x_3;x_4],[0;0]),0,[-1 1 -1 1],'b')

% store the last beta-value
bval_array(n-1)        = full(bval);

% total solver time over all iterations
solverTimes_total(n-1) = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);
buildTimes(n-1) = tempBuildTime + sum(callTime1) + sum(callTime2) + sum(callTime3);


end
end

