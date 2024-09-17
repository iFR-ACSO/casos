%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SPOTless. Implementation is based on the examples provided in the
% SPOTless toolbox and the manual
%
%--------------------------------------------------------------------------
function [bval_array, solverTimes_total,buildTimes] = roaEstNlink_benchSPOTless(Nmax,deg)




buildTimes          = zeros(Nmax-1,1);
solverTimes_total   = zeros(Nmax-1,1);
bval_array          = zeros(Nmax-1,1);

for n = Nmax

disp(['Compute maximum ROA for the ' num2str(n) '-link pendulum'])

% system states
x = msspoly ('x' , 2*n ) ;

% system dynamics
f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);

% load A and B matrix
load(['data_n' num2str(n)])


% Lyapunov function candidate
Vval = x'*P*x;
p = x'*x*10;
    
% enforce positivity
l = 1e-6*(x'*x);

% initialize arrays
endTimeBuild1 = [];
endTimeBuild2 = [];

solverTime1 = [];
solverTime2 = [];
solverTime3 = [];
bval_old = [];


% bisection tolerances 
relbistol = 1e-4;
absbistol = 1e-4;

%% V-s-iteration
for iter = 1:20

    % to make sure we do not use the old solution again
	gval = [];
    bval = [];

    % solve gamma-step

    % % find largest stable level set
    % lb = 0; ub = 1000;
    % 
    % % bisection
    % while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
    %     % trial gamma
    %     gtry = (lb+ub)/2;
        
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

         % solve problem
        opt         = spot_sdp_default_options();
        % opt.verbose = 4;
 
        sol1 = pr1.minimize(msspoly(1),@spot_mosek,opt);
      
        % sol.status
        if strcmp(string(sol1.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE") && sol1.info.mosekinfo.MSK_DINF_INTPNT_OPT_STATUS > 0.9
            % adapt lower interval bound
            % lb = gtry;

            % store latest solution
            % gval = gtry;
            s1val =  sol1.eval(s1);
        else
        %     % adapt upper interval bound
        %     ub = gtry;
            disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
            return
        end

        
        % buildTime is total time spend to setup constraints (i.e sos problem),
        % do the transcription (poly --> sdp --> poly) we subtract the
        % solver time afterwards to only consider the actual build process
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-sol1.info.wtime];
        solverTime1   = [solverTime1 sol1.info.wtime];
		
    % end


   
	% solve beta-step

    % find largest possible shape function
    lb = 0; ub = 1000;

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
    
       
        if strcmp(string(sol2.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE") && sol2.info.mosekinfo.MSK_DINF_INTPNT_OPT_STATUS > 0.9
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
	

	% solve V-step
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
	
    if strcmp(string(sol3.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE") && sol3.info.mosekinfo.MSK_DINF_INTPNT_OPT_STATUS > 0.9
           
            % extract solution
            Vval = sol3.eval(V);
	
			fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(1));
	
	
	else
		disp(['Problem infeasible in V-step in iteration:' num2str(iter)])
		break
    end


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

buildTimes(n-1)  = sum(endTimeBuild1) + sum(endTimeBuild2) + sum(endTimeBuild3);
solverTimes_total(n-1) = sum(solverTime1) + sum(solverTime2) + sum(solverTime3);
bval_array(n-1) = bval;

end % end of function


end