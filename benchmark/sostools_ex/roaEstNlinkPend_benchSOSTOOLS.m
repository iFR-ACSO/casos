%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SPOTless. Implementation is based on the examples provided in the
% SPOTless toolbox and the manual
%
%--------------------------------------------------------------------------
function [solverTimes,buildTimes,gval_array]= roaEstNlinkPend_benchSOSTOOLS(Nmax,deg)

% initialize arrays
endTimeBuild1 = [];
solverTime1 = [];

% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-3;
absbistol = 1e-3;

buildTimes  = zeros(Nmax-1,1);
solverTimes = zeros(Nmax-1,1);
gval_array  = zeros(Nmax-1,1);


%% V-s-iteration
for n = 2:Nmax
    
    % indeterminates
    x = mpvar ('x' , 2*n,1 ) ;
    
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
        prog1      = sosprogram(x);

        % decision variable
        [prog1,s1] = sossosvar(prog1,monomials(x,0:1));

        % add sos constraints
        prog1 = sosineq(prog1,s1*(V-gtry ) - jacobian(V,x)*f);

         % solve problem
        opt         = spot_sdp_default_options();
        opt.verbose = 0;
 
      
        % sol.status
        solver_opt.solver = 'mosek';
        [prog1,~] = sossolve(prog1,solver_opt);
    
        % check if feasible
        if prog1.solinfo.info.dinf == 0 && prog1.solinfo.info.pinf == 0
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
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-prog1.solinfo.info.cpusec];
        solverTime1   = [solverTime1 prog1.solinfo.info.cpusec];
        
		
    end
    gval_array(n-1)  = gval;
    buildTimes(n-1)  = sum(endTimeBuild1) ;
    solverTimes(n-1) = sum(solverTime1);


end % end for-loop



end % end of function
