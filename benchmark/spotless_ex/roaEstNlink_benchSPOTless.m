%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the N-link pendulum ROA 
% problem in SPOTless. Implementation is based on the examples provided 
% in the SPOTless toolbox and the manual
%
%--------------------------------------------------------------------------
function [bval_array, solverTimes_total,buildTimes] = roaEstNlink_benchSPOTless(Nmax,deg,noRuns)

buildTimes          = zeros(noRuns,Nmax-1,1);
solverTimes_total   = zeros(noRuns,Nmax-1,1);
bval_array          = zeros(noRuns,Nmax-1,1);


for jj = 1:noRuns
    for n = 2:Nmax
    
    disp(['Compute maximum ROA for the ' num2str(n) '-link pendulum'])

    % system states
    x = msspoly ('x' , 2*n ) ;
    
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
    endTimeBuild1 = [];
    endTimeBuild2 = [];
    
    solverTime1 = [];
    solverTime2 = [];
    solverTime3 = [];
    
    % bisection tolerances 
    relbistol = 1e-3;
    absbistol = 1e-3;
    
    %% V-s-iteration
    for iter = 1:20
        
        % to make sure we do not use the old solution again
        bval = [];
        
        %% solve s1 step
        % start time measure
        startTimeBuild1 = tic;
        
        % re-initialize a spot program
        pr1       = spotsosprog;
        pr1       = pr1.withIndeterminate(x);
    
        % decision variable
        [pr1,s1]  = pr1.newFreePoly(monomials(x,2:4));
    
        % add sos constraints
        pr1 = pr1.withSOS( s1 );
        pr1 = pr1.withSOS( s1*(Vval-1) - diff(Vval,x)*f - l );
         
        opt         = spot_sdp_default_options();
    
        % solve problem
        sol1 = pr1.minimize(msspoly(1),@spot_mosek,opt);
      
        % sol.status
        if strcmp(string(sol1.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE") 
    
            s1val =  sol1.eval(s1);
          
        else
    
            disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
            return
        end
    
            
        % buildTime is total time spend to setup constraints (i.e sos problem),
        % do the transcription (poly --> sdp --> poly) we subtract the
        % solver time afterwards to only consider the actual build process
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-sol1.info.wtime];
        solverTime1   = [solverTime1 sol1.info.wtime];
		    
	    %% solve beta-step
    
        % find largest possible shape function
        lb = -100; ub = 100;
    
        % biscetion
        while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
    
            % trial beta
            btry = (lb+ub)/2;
            
            % start time measure
            startTimeBuild2 = tic;
    
            % re-initialize sos program
            pr2       = spotsosprog;
            pr2       = pr2.withIndeterminate(x);
    
            % decision variable
            [pr2,s2]  = pr2.newFreePoly(monomials(x,0:2));
    
            % constraints    
            pr2 = pr2.withSOS( s2 );
            pr2 = pr2.withSOS( s2*(p-btry) + 1 - Vval );
            
            % solve problem
            opt         = spot_sdp_default_options();
            opt.verbose = 0;
            
            sol2 = pr2.minimize(msspoly(1),@spot_mosek,opt);
        
           
            if strcmp(string(sol2.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE") 
                % adapt lower interval bound
                lb  = btry;
        
                % store latest solution
                bval = btry;
                s2val = sol2.eval(s2);
    
          
            else
                 % adapt upper interval bound
                ub = btry;
            end
                % buildTime is total time spend to setup constraints (i.e sos problem),
                % do the transcription (poly --> sdp --> poly) we subtract the
                % solver time afterwards to only consider the actual build process
                endTimeBuild2 = [endTimeBuild2 toc(startTimeBuild2)-sol2.info.wtime];
                solverTime2   = [solverTime2 sol2.info.wtime];
        end
       
    
        if ~isempty(bval)
            % fprintf('beta is %g.\n', bval)
        else
             disp(['Problem infeasible in beta-step in iteration:' num2str(iter)])
            return
        end
	    
    
	    %% solve V-step
        % start time measure
	    startTimeBuild3 = tic;
	    
        % re-initialize sos program
        pr3       = spotsosprog;
        pr3       = pr3.withIndeterminate(x);
    
        % decision variables
        [pr3,V]   =  pr3.newFreePoly(monomials(x,2));
    
        % constraints
        pr3 = pr3.withSOS( V -l  );
        pr3 = pr3.withSOS( s1val*(V-1) - diff(V,x)*f -l );
        pr3 = pr3.withSOS( s2val*(p-bval) + 1 - V );
    
    
        % solve problem
        opt         = spot_sdp_default_options();
        opt.verbose = 0;
        
        sol3 = pr3.minimize(msspoly(1),@spot_mosek,opt);
	    
        % buildTime is total time spend to setup constraints (i.e sos problem),
        % do the transcription (poly --> sdp --> poly) we subtract the
        % solver time afterwards to only consider the actual build process
	    endTimeBuild3(iter) = toc(startTimeBuild3)-sol3.info.wtime;
	    solverTime3         = [solverTime3 sol3.info.wtime];
	    
        if strcmp(string(sol3.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE") 
               
            % extract solution
            Vval = sol3.eval(V);
    
		    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(1));
	    
	    
	    else
		    disp(['Problem infeasible in V-step in iteration:' num2str(iter)])
		    break
        end
    
    
    end % end for-loop (inner)
    
    % store the last beta-value
    bval_array(jj,n-1)        = bval;

    % total solver time over all iterations
    buildTimes(jj,n-1)        = sum(endTimeBuild1) + sum(endTimeBuild2) + sum(endTimeBuild3);
    solverTimes_total(jj,n-1) = sum(solverTime1) + sum(solverTime2) + sum(solverTime3);

    
    end % end of  for loop N-link
end % end for loop noRuns
end % end of function