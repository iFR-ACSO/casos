function buildproblem(obj,qcsos)
% Build SOS problem for bisection.

opts = obj.opts;

% detect objective
if is_symbolic(qcsos.f)
    dvar = qcsos.f;

    % positive gradient (minimization)
    obj.qc_sign = +1;
else
    dvar = -qcsos.f;

    % negative gradient (maximization)
    obj.qc_sign = -1;
end

% problem size
n = length(qcsos.x);
m = length(qcsos.g);

% get cone dimensions
Nl = get_dimension(obj.get_cones,opts.Kx,'lin');
Ns = get_dimension(obj.get_cones,opts.Kx,'sos');
Ml = get_dimension(obj.get_cones,opts.Kc,'lin');
Ms = get_dimension(obj.get_cones,opts.Kc,'sos');

assert(n == (Nl + Ns), 'Dimension of Kx must be equal to number of variables (%d).', n);
assert(m == (Ml + Ms), 'Dimension of Kc must be equal to number of constraints (%d).', m)

% select sum-of-squares variables and constraints
Is = [false(Nl,1); true(Ns,1)];
Js = [false(Ml,1); true(Ms,1)];

% build SOS problem
sos.x = qcsos.x;
sos.g = qcsos.g;
sos.f = 0;
sos.p = [qcsos.p; dvar];

% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = struct('lin',Nl,'sos',Ns);
sosopt.Kc = struct('lin',Ml,'sos',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos,sosopt);

% store basis
obj.sparsity_x  = obj.sossolver.sparsity_x;
obj.sparsity_xl = obj.sossolver.sparsity_xl;
obj.sparsity_xs = obj.sossolver.sparsity_xs;
obj.sparsity_p  = sparsity(qcsos.p);
obj.sparsity_f  = sparsity(qcsos.f);
obj.sparsity_g  = obj.sossolver.sparsity_g;
obj.sparsity_gl = obj.sossolver.sparsity_gl;
obj.sparsity_gs = obj.sossolver.sparsity_gs;

end
