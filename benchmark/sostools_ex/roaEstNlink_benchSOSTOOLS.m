%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration forthe N-link pendulum dynamics 
% ROA problem in SOSTOOLS using the default dpvar data structure. 
% Implementation is based on the examples provided in the
% SOSTOOLS toolbox and the documentation.
%
%--------------------------------------------------------------------------

function [bval_array, solverTimes_total,buildTimes] = roaEstNlink_benchSOSTOOLS(Nmax,deg,noRuns)


buildTimes          = zeros(noRuns,Nmax-1);
solverTimes_total   = zeros(noRuns,Nmax-1);
bval_array          = zeros(noRuns,Nmax-1);

for jj = 1:noRuns
    for n = 2:Nmax
        
    disp(['Compute maximum ROA for the ' num2str(n) '-link pendulum'])

    x = mpvar('x',2*n,1);
    
    % system dynamics
    f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);
    
    % load P matrix
    load(['data_n' num2str(n)])
    
    % Lyapunov function candidate
    Vval = x'*P*x;
    p = x'*x*100;
        
    % enforce positivity
    l = 1e-6*(x'*x);
    
    % initialize arrays
    bval = [];
    
    endTimeBuild1 = [];
    endTimeBuild2 = [];
    
    solverTime1 = [];
    solverTime2 = [];
    solverTime3 = [];
    
    bval_old = [];
    
    % bisection tolerances 
    relbistol = 1e-3;
    absbistol = 1e-3;
    
    %% V-s-iteration
    for iter = 1:20
    
        % to make sure we do not use the old solution again
        bval = [];
    
        %% solve s1-step
    
        % start time measure
        startTimeBuild1 = tic;
        % re-initialize sos program
        prog1      = sosprogram(x);
        
        % decision variables
        [prog1,s1] = sossosvar(prog1,monomials(x,1:2));
        
        % sos constraints
        prog1 = sosineq(prog1,s1*(Vval-1 ) - jacobian(Vval,x)*f - l);
        
        %  call solver
        solver_opt.solver = 'mosek';
        solver_opt.simplify = 1;
        [prog1,~] = sossolve(prog1,solver_opt);
    
        % check if feasible
        if prog1.solinfo.info.dinf == 0 && prog1.solinfo.info.pinf == 0 && prog1.solinfo.info.feasratio >= 0.9
            s1val = sosgetsol(prog1,s1);
        else
            disp('s1 step infeasible')
            break
        end
        
        % buildTime is total time spend to setup constraints (i.e sos problem),
        % do the transcription (poly --> sdp --> poly) we subtract the
        % solver time afterwards to only consider the actual build process
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-prog1.solinfo.info.wallTime];
        solverTime1   = [solverTime1 prog1.solinfo.info.wallTime];
	     
    
        %% solve beta-step
        lb = -100; ub = 100;
        
        % bisection
        while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
            % trial beta
            btry = (lb+ub)/2;
            
            % start time measure
            startTimeBuild2 = tic;
    
            % re-initialize sos program
            prog2      = sosprogram(x);
    
            % decision variable
            [prog2,s2] = sossosvar(prog2,monomials(x,0:1));
            
            % constraints
            prog2 = sosineq(prog2,s2*(p-btry ) + 1 - Vval);
            
            %  call solver
            solver_opt.solver   = 'mosek';
            solver_opt.simplify = 1;
            [prog2,~]           = sossolve(prog2,solver_opt);
        
      
            if prog2.solinfo.info.dinf == 0 && prog2.solinfo.info.pinf == 0 && prog2.solinfo.info.feasratio >= 0.9
                % adapt lower interval bound
                lb    = btry;
    
                % store latest solution
                bval  = btry;
                s2val = sosgetsol(prog2,s2);
            else
                % adapt upper interval bound
                ub = btry;
            end
            
            % buildTime is total time spend to setup constraints (i.e sos problem),
            % do the transcription (poly --> sdp --> poly) we subtract the
            % solver time afterwards to only consider the actual build process
            endTimeBuild2 = [endTimeBuild2 toc(startTimeBuild2)-prog2.solinfo.info.wallTime];
    
            % solver time is mosek CPU time; slightly different  to tic-toc
            % measurements
            solverTime2   = [solverTime2 prog2.solinfo.info.wallTime];
	     
        end
    
        if ~isempty(bval)
            % fprintf('gamma is %g.\n', gval)
        else
             disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
            return
        end
    
    
        %% Solve V-step
        % start time measure
	    startTimeBuild3 = tic;
	    
        % re-initialize sos program
        prog3     = sosprogram(x);
    
        % decision variable
        [prog3,V] = sospolyvar(prog3,monomials(x,2),[]);
    
        % constraints
        prog3 = sosineq(prog3,V-l);
        prog3 = sosineq(prog3,s1val*(V-1) - jacobian(V,x)*f - l);
        prog3 = sosineq(prog3,s2val*(p-bval) +  1 - V);
        
        %  call solver
        solver_opt.solver = 'mosek';
        solver_opt.simplify = 1;
        [prog3,~] = sossolve(prog3,solver_opt);
    
	    
	    endTimeBuild3(iter) = toc(startTimeBuild3)-prog3.solinfo.info.wallTime;
	    solverTime3         = [solverTime3 prog3.solinfo.info.wallTime];
	    
	    if prog3.solinfo.info.dinf == 0 && prog3.solinfo.info.pinf == 0 && prog3.solinfo.info.feasratio >= 0.9
    
		        % extract solution 
                Vval = sosgetsol(prog3,V);
	            
                % print progress
			    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(1));
			    
	    else
		    disp(['Problem infeasible in V-step in iteration:' num2str(iter)])
		    break
	    end
    
    
    end % end for-loop (inner)
    
    % store the last beta-value
    bval_array(jj,n-1)        = full(bval);
    
    % total solver time over all iterations
    solverTimes_total(jj,n-1) = sum(solverTime1) + sum(solverTime2) + sum(solverTime3);
    buildTimes(jj,n-1)        = sum(endTimeBuild1) + sum(endTimeBuild2) + sum(endTimeBuild3);
    
    
    end % end of for loop for N-link
end % end of for loop for n-runs
end % end of function