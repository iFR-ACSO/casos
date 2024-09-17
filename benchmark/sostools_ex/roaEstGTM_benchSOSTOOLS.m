%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SOSTOOLS using the default dpvar data structure. 
% Implementation is based on the examples provided in the
% SOSTOOLS toolbox and the documentation.
%
%--------------------------------------------------------------------------

function [gval,bval,solverTime,buildTime]= roaEstGTM_benchSOSTOOLS()


pvar x1 x2 x3 x4;
x= [x1; x2; x3;x4];

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

% initialize arrays
bval = [];
gval = [];

endTimeBuild1 = [];
endTimeBuild2 = [];

solverTime1 = [];
solverTime2 = [];
solverTime3 = [];

bval_old = [];

% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-3;
absbistol = 1e-3;

%% V-s-iteration
for iter = 1:100

    % to make sure we do not use the old solution again
    bval = [];
    gval = [];
    
    % solve gamma-step
    
    % find largest stable level set
    lb = 0; ub = 1000;
    
    % bisection 
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
            
        % trial gamma
        gtry = (lb+ub)/2;

        % start time measure
        startTimeBuild1 = tic;

        % re-initialize sos program
        prog1      = sosprogram(x);
        
        % decision variables
        [prog1,s1] = sossosvar(prog1,monomials(x,1:2));
        
        % sos constraints
        prog1 = sosineq(prog1,s1*(Vval-gtry ) - jacobian(Vval,x)*f - l);
        
        %  call solver
        solver_opt.solver = 'mosek';
        solver_opt.simplify = 1;
        [prog1,~] = sossolve(prog1,solver_opt);
    
        % check if feasible
        if prog1.solinfo.info.dinf == 0 && prog1.solinfo.info.pinf == 0

            % adapt lower interval bound
            lb   = gtry;

            % store latest solution
            gval  = gtry;
            s1val = sosgetsol(prog1,s1);
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



    if ~isempty(gval)
        % fprintf('gamma is %g.\n', gval)
    else
         disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
        return
    end


    % solve beta-step
    lb = 0; ub = 1000;
    
    % bisection
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
        % trial beta
        btry = (lb+ub)/2;
        
        % start time measure
        startTimeBuild2 = tic;

        % re-initialize sos program
        prog2      = sosprogram(x);

        % decision variable
        [prog2,s2] = sossosvar(prog2,monomials(x,0:2));
        
        % constraints
        prog2 = sosineq(prog2,s2*(p-btry ) + gval - Vval);
        
        %  call solver
        solver_opt.solver = 'mosek';
        solver_opt.simplify = 1;
        [prog2,~]     = sossolve(prog2,solver_opt);
    
    
        if prog2.solinfo.info.dinf == 0 && prog2.solinfo.info.pinf == 0
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
        endTimeBuild2 = [endTimeBuild1 toc(startTimeBuild2)-prog2.solinfo.info.wallTime];

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


    %Solve V-step
    % start time measure
	startTimeBuild3 = tic;
	
    % re-initialize sos program
    prog3     = sosprogram(x);

    % decision variable
    [prog3,V] = sospolyvar(prog3,monomials(x,2:4),[]);

    % constraints
    prog3 = sosineq(prog3,V-l);
    prog3 = sosineq(prog3,s1val*(V-gval) - jacobian(V,x)*f - l);
    prog3 = sosineq(prog3,s2val*(p-bval) +  gval - V);
    
    %  call solver
    solver_opt.solver = 'mosek';
    solver_opt.simplify = 1;
    [prog3,~] = sossolve(prog3,solver_opt);

	
	endTimeBuild3(iter) = toc(startTimeBuild3)-prog3.solinfo.info.wallTime;
	solverTime3         = [solverTime3 prog3.solinfo.info.wallTime];
	
	if prog3.solinfo.info.dinf == 0 && prog3.solinfo.info.pinf == 0

		    % extract solution 
            Vval = sosgetsol(prog3,V);
	        
            % print progress
			fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(gval));
			
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

buildTime  = sum(endTimeBuild1) + sum(endTimeBuild2) + sum(endTimeBuild3);
solverTime = sum(solverTime1)   + sum(solverTime2)   + sum(solverTime3);
 

% save the complete workspace, so people do not have to re-run execpt they
% want to
% save('SOSTOOLS_GTM_ROA_bench.mat')

end