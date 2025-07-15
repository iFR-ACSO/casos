%--------------------------------------------------------------------------
% 
% Implementation of sequential SOSS for the N-link pendulum dynamics
% in CaSoS.
%
%--------------------------------------------------------------------------

function [fval_array, solverTimes_total,buildTimes,iteration_array,solStatus_array,solverStats] = bench_Nlink_Seq_badInit(noRuns,deg,Nmax)

import casos.toolboxes.sosopt.*

buildTimes          = zeros(noRuns,Nmax-1);
solverTimes_total   = zeros(noRuns,Nmax-1);
fval_array          = zeros(noRuns,Nmax-1);



% compute each N-link noRuns times
for jj = 1:noRuns

    % iterate of the Number of links
    for n = 2:Nmax
    
    disp(['Compute maximum ROA for the ' num2str(n) '-link pendulum'])
    
    x = casos.PS('x',2*n,1);
    
    % get dynamics
    f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);
    
    % load P matrix
    load(['data_n' num2str(n)],'P');
    
    % Lyapunov function candidate
    Vval = x'*P*x;
 
    % enforce positivity
    l = 1e-6*(x'*x);
    
    % start to measure build time of all parameterized solver
    buildTimes_start = tic;
    
    % Lyapunov function candidate
    V = casos.PS.sym('v',monomials(x,2));
    
    % SOS multiplier
    s1 = casos.PS.sym('s1',monomials(x,2));
    
    
    
    g = x'*P*0.01*x-2;

    % fixed gamma
    gval = 1;
    
    % % 
    % figure()
    % pcontour(subs(g, x(3:end), zeros( length(x(3:end)),1)),0,[-1 1 -1 1],'r')
    % hold on
    % pcontour(subs(Vval, x(3:end), zeros( length(x(3:end)),1)),1,[-1 1 -1 1])
    % pcontour(Vval,1,[-6 6 -6 6])
    %% setup solve
    % solver 3: V-step
    sos = struct('x',[V;s1], ...    % dec.var
                  'f',dot(g-(V),g-(V)));
    
    % constraints
    sos.('g') = [V-l; 
                  s1;
                  s1*(V-gval)-nabla(V,x)*f-l];
    
    opts     = struct;
    opts.Kx = struct('lin', 2); 
    opts.Kc = struct('sos', 3);
    opts.verbose = 1;
    % opts.tolerance_opt = 1e-3;
    opts.sossol                     = 'mosek';

    % build third solver
    S = casos.nlsossol('S','sequential',sos,opts);
    
    tempBuildTime = toc(buildTimes_start);

    
    sol = S('x0',casos.PD([-(x'*x),-(x'*x)]));
    

    % store the last cost-function-value
    fval_array(jj,n-1)        = full(sol.f);
    iteration_array(jj,n-1)   = length(S.stats.single_iterations);
    solStatus_array{jj,n-1}   = S.stats.solutionStatus;
    solverStats{jj,n-1}       = {S.stats.single_iterations};
    % total solver time over all iterations
    solverTimes_total(jj,n-1) = S.stats.totalSolveTime;
    buildTimes(jj,n-1)        = tempBuildTime;
     
    end % end for-loop (inner)
    
    
end % end of for loop for noRuns
end % end-of-function