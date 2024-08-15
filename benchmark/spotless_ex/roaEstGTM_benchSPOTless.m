function [gval,bval,solveTime,solverTime,buildTime]= roaEstGTM_benchSPOless()

x = msspoly ('x' , 4 ) ;

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
        
        pr1       = spotsosprog;
        pr1       = pr1.withIndeterminate(x);
        [pr1,s1]  = pr1.newFreePoly(monomials(x,2:4));


        pr1 = pr1.withSOS( s1 );
        pr1 = pr1.withSOS( s1*(Vval-gtry) - diff(Vval,x)*f - l );

        opt = spot_sdp_default_options();
        opt.verbose = 0;
        

        sol1 = pr1.minimize(msspoly(1),@spot_mosek,opt);
      
        % sol.status
        if strcmp(string(sol1.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE")
            lb = gtry;
            gval = gtry;
            s1val =  sol1.eval(s1);
        else
            ub = gtry;
        end

        
        % buildtime is complete time - solver time 
        endTimeParse1 = [endTimeParse1 toc(startTimeParse1)-sol1.info.wtime];
        solverTime1   = [solverTime1 sol1.info.wtime];
		
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
    
        pr2       = spotsosprog;
        pr2       = pr2.withIndeterminate(x);
        [pr2,s2]  = pr2.newFreePoly(monomials(x,0:4));


        pr2 = pr2.withSOS( s2 );
        pr2 = pr2.withSOS( s2*(p-btry) + gval - Vval );

        opt = spot_sdp_default_options();
        opt.verbose = 0;
        

        sol2 = pr2.minimize(msspoly(1),@spot_mosek,opt);
    
       
        if strcmp(string(sol2.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE")
            lb  = btry;
    
            % make sure to keep the feasible solution
            bval = btry;
            s2val = sol2.eval(s2);
            
        else
            ub = btry;
        end


            endTimeParse2 = [endTimeParse2 toc(startTimeParse2)-sol2.info.wtime];
            solverTime2   = [solverTime2 sol2.info.wtime];
    end
   

    if ~isempty(bval)
        % fprintf('beta is %g.\n', bval)
    else
         disp(['Problem infeasible in beta-step in iteration:' num2str(iter)])
        return
    end
	

	% solve V-step
	startTimeParse3 = tic;
	
        pr3       = spotsosprog;
        pr3       = pr3.withIndeterminate(x);
        [pr3,V]  =  pr3.newFreePoly(monomials(x,2:4));


        pr3 = pr3.withSOS( V -l  );
        pr3 = pr3.withSOS( s2val*(p-bval) + gval - V );
        pr3 = pr3.withSOS( s1val*(V-gval) -diff(V,x)*f -l );

        opt = spot_sdp_default_options();
        opt.verbose = 0;
        

        sol3 = pr3.minimize(msspoly(1),@spot_mosek,opt);
	
	endTimeParse3(iter) = toc(startTimeParse3)-sol3.info.wtime;
	solverTime3         = [solverTime3 sol3.info.wtime];
	
	if strcmp(string(sol3.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE")
				Vval = sol3.eval(V);
	
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
