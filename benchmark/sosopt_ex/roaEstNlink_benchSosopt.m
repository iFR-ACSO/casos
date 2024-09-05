%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SOSOPT/GSOSOPT. GSOSOPT performs the bisection.
% Implementation is based on the example gsosoptdemo1.m provide in the
% SOSOPT toolbox.
%
%--------------------------------------------------------------------------
function [bval_array, solverTimes_total,buildTimes] = roaEstNlink_benchSosopt(Nmax,deg)

buildTimes          = zeros(Nmax-1,1);
solverTimes_total   = zeros(Nmax-1,1);
bval_array          = zeros(Nmax-1,1);

for n = 2:Nmax
x = mpvar('x',2*n,1);

% system dynamics
f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);

% load A and B matrix
load(['data_n' num2str(n)])


% Lyapunov function candidate
Vval = x'*P*x;
p = x'*x*100;

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
V  = polydecvar('v',monomials(x, 2 ));
s1 = sosdecvar('s1',monomials(x, 1:2 ));
s2 = sosdecvar('s2',monomials(x, 0:2 ));

%% V-s-iteration
for iter = 1:20
  
    opts = sosoptions;
    opts.solver = 'mosek';

    % gamma-step    
    startTimeBuild1 = tic;

    sosc    = polyconstr;
    sosc(1) = s1 >=0;
    % different sign compared to paper; just another definition in
    % bisection
    sosc(2) = s1*(Vval -1) - jacobian(Vval,x)*f - l >= 0;
    
    % Solve with gsosopt
    [info,dopt] = sosopt(sosc,x,opts);

     % extract solution
    s1val = subs(s1,dopt);

    % time measure of first subproblem
    solverTime1(iter) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); in gsosopt the problem
    % is automatically re-build in each subiteration; we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild1     = [endTimeBuild1 toc(startTimeBuild1)-solverTime1(iter)]; 


    % beta-step    
     opts = gsosoptions;
    opts.solver = 'mosek';

    opts.minobj = -1000; 
    opts.maxobj = 0;
    startTimeBuild2 = tic;
    sosc    = polyconstr;

    sosc(1) = s2 >=0;
    
  
    % different sign compared to paper; just another definition in
    % bisection
    sosc(2) = s2*(p + b) + 1 - Vval >= 0;
    
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

    % beta-step    
    startTimeBuild3 = tic;

    sosc    = polyconstr;
    sosc(1) = V - l >=0;
    sosc(2) = s1val*(V - 1) - jacobian(V,x)*f - l >= 0;
    sosc(3) = s2val*(p - bval) + 1 - V >= 0;
    
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

    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(1));
	
				
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

end % end-for-loop

% store the last beta-value
bval_array(n-1)        = full(bval);

% total solver time over all iterations
solverTimes_total(n-1) = sum(solverTime1)   + sum(solverTime2)   + sum(solverTime3);
buildTimes(n-1)        = sum(endTimeBuild1) + sum(endTimeBuild2) + sum(endTimeBuild3);

end


end % end of function