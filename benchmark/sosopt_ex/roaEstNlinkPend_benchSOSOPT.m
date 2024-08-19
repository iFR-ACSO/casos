%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SPOTless. Implementation is based on the examples provided in the
% SPOTless toolbox and the manual
%
%--------------------------------------------------------------------------
function [solverTimes,buildTimes,gval_array]= roaEstNlinkPend_benchSOSOPT(Nmax,deg)


% initialize arrays
endTimeBuild1 = [];
solverTime1 = [];
gval_array  = zeros(Nmax-1,1);
pvar g;

%% V-s-iteration
for n = 2:Nmax
    
    % indeterminates
    x  = mpvar('x',[2*n,1]);
    s = sosdecvar('s',monomials(x, 0:1 ));
    
    % system dynamics
    f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);
    load(['data_n' num2str(n)])
    % Lyapunov function candidate
    V = x'*S*x;
    p = x'*x;

    opts = gsosoptions;
    opts.minobj = -100; 
    opts.maxobj = 100;
    opts.solver = 'mosek';

    % gamma-step    
    startTimeBuild1 = tic;

    sosc    = polyconstr;
    sosc(1) = s >=0;
    % different sign compared to paper; just another definition in
    % bisection
    sosc(2) = s*(V + g) - jacobian(V,x)*f >= 0;
    
    % Solve with gsosopt
    [info,dopt] = gsosopt(sosc,x,g,opts);

     % extract solution
    % s1val = subs(s1,dopt);
    gval  = -info.tbnds(2); 
    
    % time measure of first subproblem
    solverTimes(n-1) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); in gsosopt the problem
    % is automatically re-build in each subiteration; we subtract the
    % solver time afterwards to only consider the actual build process
     buildTimes(n-1)     = [endTimeBuild1 toc(startTimeBuild1)-solverTimes(n-1)]; 
     gval_array(n-1)   = gval;
end

end % end of function

