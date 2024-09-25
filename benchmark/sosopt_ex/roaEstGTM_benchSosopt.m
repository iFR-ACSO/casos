%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SOSOPT/GSOSOPT. GSOSOPT performs the bisection.
% Implementation is based on the example gsosoptdemo1.m provide in the
% SOSOPT toolbox.
%
%--------------------------------------------------------------------------
function [gval,bval,solverTime,buildTime]= roaEstGTM_benchSosopt()

% System dynamics
pvar x1 x2 x3 x4

x = [x1; x2; x3;x4];

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
endTimeBuild1 = [];
endTimeBuild2 = [];
bval_old = [];

solverTime1 = zeros(100,1);
solverTime2 = zeros(100,1);
solverTime3 = zeros(100,1);

% decision variables
pvar g b;
V  = polydecvar('v',monomials(x, 2:4 ));
s1 = sosdecvar('s1',monomials(x, 1:2 ));
s2 = sosdecvar('s2',monomials(x, 0:2 ));

%% V-s-iteration
for iter = 1:100
  
    opts = gsosoptions;
    opts.minobj = -1000; 
    opts.maxobj = 0;
    opts.solver = 'mosek';

    %% gamma-step    
    startTimeBuild1 = tic;

    sosc    = polyconstr;
    sosc(1) = s1 >=0;
    % different sign compared to paper; just another definition in
    % bisection
    sosc(2) = s1*(Vval + g) - jacobian(Vval,x)*f - l >= 0;
    
    % Solve with gsosopt
    [info,dopt] = gsosopt(sosc,x,g,opts);

     % extract solution
    s1val = subs(s1,dopt);
    gval  = -info.tbnds(2); 
    
    % time measure of first subproblem
    solverTime1(iter) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); in gsosopt the problem
    % is automatically re-build in each subiteration; we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild1     = [endTimeBuild1 toc(startTimeBuild1)-solverTime1(iter)]; 


    %% beta-step    
    startTimeBuild2 = tic;
    sosc    = polyconstr;
    sosc(1) = s2 >=0;

    % different sign compared to paper; just another definition in
    % bisection
    sosc(2) = s2*(p + b) + gval - Vval >= 0;
    
    % Solve with gsosopt
    [info,dopt] = gsosopt(sosc,x,b,opts);
    
     % extract solution
    s2val = subs(s2,dopt);
    bval  = -info.tbnds(2); 
    
    % time measure of second subproblem
    solverTime2(iter) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); in gsosopt the problem
    % is automatically re-build in each subiteration; we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild2    = [endTimeBuild2 toc(startTimeBuild2)-solverTime2(iter)];

    %% V-step    
    startTimeBuild3 = tic;

    sosc    = polyconstr;
    sosc(1) = V - l >=0;
    sosc(2) = s1val*(V - gval) - jacobian(V,x)*f - l >= 0;
    sosc(3) = s2val*(p - bval) + gval - V >= 0;
    
    % Solve with sosopt
    opts = sosoptions;
    opts.solver = 'mosek';
      
    [info,dopt,~] = sosopt(sosc,x,opts);

    % time measure of third subproblem
    solverTime3(iter)   = info.sdpsol.solverinfo.solvertime;

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild3(iter) = toc(startTimeBuild3)-solverTime3(iter);
    
    % extract solution
    Vval = subs(V,dopt);

    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(gval));
	
				
	%% check convergence
	if ~isempty(bval_old)
			if abs(full(bval-bval_old)) <= 1e-3
				break
			else
				bval_old = bval;
			end
		else
			bval_old = bval;
	end

end % end-for-loop

%% prepare output
buildTime  = sum(endTimeBuild1) + sum(endTimeBuild2) + sum(endTimeBuild3);
solverTime = sum(solverTime1)   + sum(solverTime2)   + sum(solverTime3);


% save the complete workspace, so people do not have to re-run execpt they
% want to
% save('SOSOPT_GTM_ROA_bench.mat')

end % end of function