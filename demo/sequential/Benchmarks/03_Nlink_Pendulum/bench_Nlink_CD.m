%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the N-link pendulum dynamics
% in CaSoS.
%
%--------------------------------------------------------------------------

function [fval_array, solverTimes_total,buildTimes,iteration_array,solverStats] = bench_Nlink_CD(noRuns,deg,Nmax)


buildTimes          = zeros(noRuns,Nmax-1);
solverTimes_total   = zeros(noRuns,Nmax-1);
bval_array          = zeros(noRuns,Nmax-1);



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
    
    
    % build a safe set
    g = x'*P*0.01*x-2;

    % fixed gamma
    gval = 1;
    
    %% setup solver
    
    % solver 1: s1-step
    sos1 = struct('x',s1,'p',V);     % parameter
    
    % constraint
    sos1.('g') = [s1; s1*(V-gval)-nabla(V,x)*f-l];
    
    % states + constraint are SOS cones
    opts.Kx = struct('lin', 1);
    opts.Kc = struct('sos', 2);
    
    % build first solver
    S1 = casos.sossol('S1','mosek',sos1,opts);
    
    % solver 3: V-step
    sos2 = struct('x',V, ...        % dec.var
                  'f',dot(g-V,g-V),...
                  'p',s1); % parameter
    
    % constraints
    sos2.('g') = [V-l; 
                  s1*(V-gval)-nabla(V,x)*f-l];
    
    opts     = struct;
    opts.Kx = struct('sos', 0, 'lin', 1); 
    opts.Kc = struct('sos', 2);
    
    % build third solver
    S2 = casos.sossol('S','mosek',sos2,opts);
    
    tempBuildTime = toc(buildTimes_start);
    
    % initialize arrays
    solvetime_all1 = zeros(100,1);
    solvetime_all2 = zeros(100,1);
    solvetime_all3 = zeros(100,1);
    
    % callTime1 = zeros(100,1);
    % callTime2 = zeros(100,1);
    % callTime3 = zeros(100,1);
    
    fval_old = [];
    % 
    % casos.Function('f',{poly2basis(x_full)})
    % Vval = -(x'*x);
    %% V-s-iteration
    for iter = 1:100
    
        %% s1-step
        s1start = tic;
         
        % call solver
        sol1 = S1('p',Vval);    
    
        solvetime_all1(iter) =  toc(s1start); S1.stats.mosek_info.MSK_DINF_OPTIMIZER_TIME;
        % callTime1(iter)      = toc(s1start)-solvetime_all1(iter);
    
        % extract solution
        s1val = sol1.x;
    
        %% V-step
        S2start = tic;
    
        sol2 = S2('p',s1val);
    
        solvetime_all3(iter) = toc(S2start);%S2.stats.mosek_info.MSK_DINF_OPTIMIZER_TIME;
        % callTime3(iter) = toc(S2start)-solvetime_all3(iter);
    
        % extract solution
        Vval = sol2.x;
    
        % show progress 
        fprintf('Iteration %d: f = %g, g = %g.\n',iter,full(sol2.f),full(1));
        

     if ~isempty(fval_old)
        if abs(full(sol2.f-fval_old)) <= 1e-3
            break
        else
            fval_old = sol2.f;
        end
    else
        fval_old = sol2.f;
    end


        
    end % end for-loop (inner)
     
    % store the last beta-value
    fval_array(jj,n-1)        = full(sol2.f);
    iteration_array(jj,n-1)   = iter;
    % solStatus_array{jj,n-1}   = S2.stats.solutionStatus;
    solverStats{jj,n-1}       = {{S1.stats.conic,S2.stats.conic}};
    % total solver time over all iterations
    solverTimes_total(jj,n-1) = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);
    buildTimes(jj,n-1)        = tempBuildTime;% + sum(callTime1) + sum(callTime2) + sum(callTime3);
    
    
    end % end of for loop for N-link
end % end of for loop for noRuns
end % end-of-function