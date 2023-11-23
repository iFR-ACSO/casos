% Validate stability with LF candidate.

% system states
x = casos.PS('x',2,1);
% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];
% Lyapunov function candidate
V = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;
% derivative w.r.t. time
Vdot = jacobian(V,x)*f;
% SOS multiplier
s = casos.PS.sym('q',monomials(x,0:1),'gram');
% enforce positivity
l = 1e-6*(x'*x);
% level of stability
g = casos.PS.sym('g');

% define SOS feasibility
sos = struct('x',s,'g',[V-l; s*(V-g)-Vdot-l],'p',g);
% states + constraint are SOS cones
opts.Kx.s = 1; opts.Kc.s = 2;
% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);

% find largest stable level set
lb = 0; ub = 10;
% bisection
while (ub - lb > 1e-1)
    ptry = (lb+ub)/2;

    % evaluate parametrized SOS problem
    sol = S('p',ptry);

    switch (S.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS', lb = ptry;
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}, ub = ptry;
        otherwise, error('Failed.')
    end
end

fprintf('Maximal stable level set is %g.\n', lb);
