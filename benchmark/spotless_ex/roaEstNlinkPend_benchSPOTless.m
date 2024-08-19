%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SPOTless. Implementation is based on the examples provided in the
% SPOTless toolbox and the manual
%
%--------------------------------------------------------------------------
function [solverTimes,buildTimes,gval_array]= roaEstNlinkPend_benchSPOTless(Nmax,deg)


% initialize arrays
endTimeBuild1 = [];
solverTime1 = [];


% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-3;
absbistol = 1e-3;
gval_array  = zeros(Nmax-1,1);
%% V-s-iteration
for n = 2:Nmax
    
    % indeterminates
    x = msspoly ('x' , 2*n ) ;
    
    % system dynamics
    f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);
    load(['data_n' num2str(n)])
    % Lyapunov function candidate
    V = x'*S*x;
    p = x'*x;

    % find largest stable level set
    lb = -100; ub = 100;
    
    % bisection
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
        % trial gamma
        gtry = (lb+ub)/2;
        
        % start time measure
        startTimeBuild1 = tic;
        
        % re-initialize a spot program
        pr1       = spotsosprog;
        pr1       = pr1.withIndeterminate(x);

        % decision variable
        [pr1,s]  = pr1.newFreePoly(monomials(x,0:2));

        % add sos constraints
        pr1 = pr1.withSOS( s );
        pr1 = pr1.withSOS( s*(V-gtry)-diff(V,x)*f );

         % solve problem
        opt         = spot_sdp_default_options();
        opt.verbose = 0;
 
        sol1 = pr1.minimize(msspoly(1),@spot_mosek,opt);
      
        % sol.status
        if strcmp(string(sol1.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE")
            % adapt lower interval bound
            lb = gtry;

            % store latest solution
            gval = gtry;
            % s1val =  sol1.eval(s1);
        else
            % adapt upper interval bound
            ub = gtry;
        end

        
        % buildTime is total time spend to setup constraints (i.e sos problem),
        % do the transcription (poly --> sdp --> poly) we subtract the
        % solver time afterwards to only consider the actual build process
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-sol1.info.wtime];
        solverTime1   = [solverTime1 sol1.info.wtime];
        
		
    end
    gval_array(n-1)  = gval;
    buildTimes(n-1)  = sum(endTimeBuild1) ;
    solverTimes(n-1) = sum(solverTime1);


end % end for-loop



end % end of function
