function [gval,bval,solveTime,solverTime,buildTime]= roaEstGTM_benchYALMIP(defaultOpts)

sdpvar x1 x2 x3 x4 g b

x = [x1;x2;x3;x4];

% Polynomial Dynamics
f = GTM_dynamics(x(1),x(2),x(3),x(4));

% shape function
p = x'*x*1e2;

% initial Lyapunov functionn
P = [395.382671347059	-23.0032507976836	3.16965275615691	29.2992065909380
    -23.0032507976836	90.4764915483638	-16.1191428789579	-132.594376986429
    3.16965275615691	-16.1191428789579	3.44002648214490	24.2292666551058
    29.2992065909380	-132.594376986429	24.2292666551058	202.114797577027];

Vval = x'*P*x;


l = 1e-6*(x'*x);

% polynomial(indet, maxdeg, mindeg)
[V,cv1]  = polynomial(x,4,2);
[s1,c1]  = polynomial(x,4,2);
[s2,c2]  = polynomial(x,2,0);

% use default options
% if defaultOpts
solverset = sdpsettings('solver','mosek', ...
                        'verbose',0);        
% else
% solverset = sdpsettings('solver','mosek', ...
%                         'verbose',0, ...
%                          'sos.traceobj',0,...   % Minimize trace of Gram matrix in problems without objective function
%                          'sos.newton',0,...     % Use Newton polytope to reduce size
%                          'sos.congruence',2,... % Block-diagonalize using congruence classes
%                          'sos.scale',0);        % scale polynomials
% % end

%--------------------------------------------------------------------------
% see https://yalmip.github.io/tutorial/sumofsquaresprogramming/ on how to
% setup constraint sos problems
%--------------------------------------------------------------------------


endTimeParse1 = [];
endTimeParse2 = [];

solverTime1 = [];
solverTime2 = [];
solverTime3 = [];
bval_old = [];

solveTime_start = tic;
for iter = 1:100
    			% to make sure we do not use the old solution again
				gval = [];
				bval = [];

    % solve gamma-step

    % find largest stable level set
    lb = 0; ub = 1000;
    

    % bisection
    % checkStart = tic;
    relbistol = 1e-3;
    absbistol = 1e-3;
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
        startTimeParse1 = tic;
        gtry = (lb+ub)/2;
        
        
        con1 = [sos(s1*(Vval-gtry) - jacobian(Vval,x)*f - l)
                sos(s1)];

        
        sol1 = solvesos(con1,[],solverset,c1 );

        if sol1.problem == 0
            lb = gtry;
            gval = gtry;
            s1val = replace(s1,c1,value(c1));
        else
            ub = gtry;
        end

        
        % buildtime is complete time - solver time 
        endTimeParse1 = [endTimeParse1 toc(startTimeParse1)-sol1.solvertime];
        solverTime1   = [solverTime1 sol1.solvertime];
		
    end
    % toc(checkStart)


    if ~isempty(gval)
        % fprintf('gamma is %g.\n', gval)
    else
         disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
        return
    end
   
    

	% solve beta-step

    % find largest possible shape function
    lb = 0; ub = 1000;

    relbistol = 1e-3;
    absbistol = 1e-3;
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
        startTimeParse2 = tic;
        btry = (lb+ub)/2;
    
        % parse problem
        con2 = [sos(s2*(p-btry) +  gval - Vval)
                sos(s2)];
       
        % solve problem
        sol2 = solvesos(con2,[],solverset,c2 );
    
       
        if sol2.problem == 0
            lb  = btry;
    
            % make sure to keep the feasible solution
            bval = btry;
            s2val = replace(s2,c2,value(c2));
            
        else
            ub = btry;
        end


            endTimeParse2 = [endTimeParse2 toc(startTimeParse2)-sol2.solvertime];
            solverTime2   = [solverTime2 sol2.solvertime];
    end
   

    if ~isempty(bval)
        % fprintf('beta is %g.\n', bval)
    else
         disp(['Problem infeasible in beta-step in iteration:' num2str(iter)])
        return
    end
	

	% solve V-step
	startTimeParse3 = tic;
	
	con3 = [sos(V-l)
			sos(s1val*(V-gval) - jacobian(V,x)*f - l)
			sos(s2val*(p-bval) +  gval - V)
			];
	
	
	% solve for Lyapunov function
	sol3 = solvesos(con3,[],solverset,cv1);
	
	endTimeParse3(iter) = toc(startTimeParse3)-sol3.solvertime;
	solverTime3         = [solverTime3 sol3.solvertime];
	
	if sol3.problem == 0
				Vval = replace(V,cv1,double(cv1));
	
			fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(gval));
	
	
	else
		disp(['Problem infeasible in V-step in iteration:' num2str(iter)])
		break
	end

   if ~isempty(bval_old)
        if abs(full(bval-bval_old)) <= 1e-3
            break
        else
            bval_old = bval;
        end
    else
        bval_old = bval;
    end


end

solveTime  = toc(solveTime_start);
buildTime  = sum(endTimeParse1) + sum(endTimeParse2) + sum(endTimeParse3);
solverTime = sum(solverTime1) + sum(solverTime2) + sum(solverTime3);
 

end
