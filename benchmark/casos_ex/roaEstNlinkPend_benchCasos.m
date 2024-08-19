%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% CaSoS. The quasi-convex solver is used to perform the bisections.
%
%--------------------------------------------------------------------------

function [solverTimes,buildTimes,gval_array] = roaEstNlinkPend_benchCasos(Nmax,deg)


buildTimes  = zeros(Nmax-1,1);
solverTimes = zeros(Nmax-1,1);
gval_array  = zeros(Nmax-1,1);

for n = 2:Nmax
    % system states
    x = casos.PS('x',2*n,1);
    
    % system dynamics
    f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);
    load(['data_n' num2str(n)])

    % Lyapunov function candidate
    V = x'*S*x;
    p = x'*x;
    
    buildTimeStart = tic;
    % scalar decision variable
    g = casos.PS.sym('g');
    
    % multiplier (idk)
    s = casos.PS.sym('s',monomials(x,0:1),'gram');
    
    % define SOS problem:
    sos = struct('x',s,'f',-g,'g', s*(V-g)-nabla(V,x)*f);
    
    % options
    opts = struct('sossol','mosek');
    opts.tolerance_abs = 1e-3;
    opts.tolerance_rel = 1e-3;
    opts.conf_interval = [-1000 0];
    
    % constraint is scalar SOS cone
    opts.Kx = struct('sos', 1);
    opts.Kc = struct('sos', 1);
    
    % solve by relaxation to SDP
    S1 = casos.qcsossol('S1','bisection',sos,opts);
    buildTimes(n-1) = toc(buildTimeStart);

    % evaluate
    sol = S1();
    solverTimes(n-1) = sum(cellfun(@(x) x.solvetime_matlab, S1.stats.iter));
    gval_array(n-1)  = full(-sol.f);
end

end % end of function