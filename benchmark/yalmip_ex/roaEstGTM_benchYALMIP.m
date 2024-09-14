%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% YALMIP. The implementation in Yalmip follow the following YALMIP 
% tutorial:
% https://yalmip.github.io/tutorial/sumofsquaresprogramming/ on how to
% setup constraint sos problems. 
% More specifically the section "Constrained polynomial optimization"
%
%--------------------------------------------------------------------------

function [gval,bval,solverTime,buildTime]= roaEstGTM_benchYALMIP()


% indeterminates 
sdpvar x1 x2 x3 x4 g b

x = [x1;x2;x3;x4];

% it takes very long for sdpvar to build up dyanmics!
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


% polynomial(indet, maxdeg, mindeg)
[V,cv1]  = polynomial(x,4,2);
[s1,c1]  = polynomial(x,4,2);
[s2,c2]  = polynomial(x,4,0);

% use default options except from verbosity and select mosek as solver
solverset = sdpsettings('solver','mosek', ...
                        'verbose',0);        


% solverset = sdpsettings('solver','mosek', ...
%                         'sos.newton',0,...
%                         'sos.congruence',0,...
%                         'sos.scale',0,...
%                         'verbose',0);        


% setup arrays
endTimeBuild1 = [];
endTimeBuild2 = [];


solverTime1 = [];
solverTime2 = [];
solverTime3 = [];
bval_old    = [];

% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-3;
absbistol = 1e-3;


%% V-s-iteration
for iter = 1:100
    
    % to make sure we do not use the old solution again
	gval = [];
    bval = [];

    % solve gamma-step
    lb = 0; ub = 1000;
    
    % bisection
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))

        % trial gamma
        gtry = (lb+ub)/2;

        % start time measure
        startTimeBuild1 = tic;
    
        % re-initialize sos program
        con1 = [sos(s1*(Vval-gtry) - jacobian(Vval,x)*f - l);
                sos(s1)];

        %  solve problem
        sol1 = solvesos(con1,[],solverset,c1);

        if sol1.problem == 0

            % adapt lower interval bound
            lb = gtry;

            % store latest solution
            gval  = gtry;
            s1val = replace(s1,c1,value(c1));

        else
            % adapt upper interval bound
            ub = gtry;
        end


         % buildTime is total time spend to setup constraints (i.e sos problem),
         % do the transcription (poly --> sdp --> poly) we subtract the
         % solver time afterwards to only consider the actual build process
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-sol1.solvertime];
        solverTime1   = [solverTime1 sol1.solvertime];
		
    end


    if ~isempty(gval)
        % fprintf('gamma is %g.\n', gval)
    else
         disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
        return
    end
   
    

	% solve beta-step

    % find largest possible shape function
    lb = 0; ub = 1000;

    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
        
        % trial beta
        btry = (lb+ub)/2;
        
        % start time measure
        startTimeBuild2 = tic;
    
         % re-initialize sos program
        con2 = [sos(s2*(p-btry) +  gval - Vval);
                sos(s2)];
       
        % solve problem
        sol2 = solvesos(con2,[],solverset,c2 );
    
       
        if sol2.problem == 0
            % adapt lower interval bound
            lb  = btry;
    
            % store latest solution
            bval = btry;
            s2val = replace(s2,c2,value(c2));
            
        else
            % adapt upper interval bound
            ub = btry;
        end

            % buildTime is total time spend to setup constraints (i.e sos problem),
            % do the transcription (poly --> sdp --> poly) we subtract the
            % solver time afterwards to only consider the actual build process
            endTimeBuild2 = [endTimeBuild2 toc(startTimeBuild2)-sol2.solvertime];
            solverTime2   = [solverTime2 sol2.solvertime];
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
	con3 = [sos(V-l);
			sos(s1val*(V-gval) - jacobian(V,x)*f - l);
			sos(s2val*(p-bval) +  gval - V)
			];
	
	
	% solve sos problem
	sol3 = solvesos(con3,[],solverset,cv1);
	
	endTimeBuild3(iter) = toc(startTimeBuild3)-sol3.solvertime;
	solverTime3         = [solverTime3 sol3.solvertime];
	
	if sol3.problem == 0
            % extract solution 
		    Vval = replace(V,cv1,double(cv1));

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
% save('YALMIP_GTM_ROA_bench.mat')

end % end of function
