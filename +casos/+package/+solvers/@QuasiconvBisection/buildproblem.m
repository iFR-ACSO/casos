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
% get cone dimensions for the decision variables
Nl = get_dimension(obj.get_cones,opts.Kx,'lin');
Ns = get_dimension(obj.get_cones,opts.Kx,'sos');
Nds  = get_dimension(obj.get_cones,opts.Kx,'dsos');
Nsds = get_dimension(obj.get_cones,opts.Kx,'sdsos');

% get cone dimensions for the constraints
Ml = get_dimension(obj.get_cones,opts.Kc,'lin');
Ms = get_dimension(obj.get_cones,opts.Kc,'sos');
Mds  = get_dimension(obj.get_cones,opts.Kc,'dsos');
Msds = get_dimension(obj.get_cones,opts.Kc,'sdsos');

assert(n == (Nl + Ns + Nds + Nsds), 'Dimension of Kx must be equal to number of variables (%d).', n);
assert(m == (Ml + Ms + Mds + Msds), 'Dimension of Kc must be equal to number of constraints (%d).', m)

% select sum-of-squares variables and constraints
Is = [false(Nl,1); true(Ns,1); true(Nds,1); true(Nsds,1)];
Js = [false(Ml,1); true(Ms,1); true(Mds,1); true(Msds,1)];

% build SOS problem
sos.x = qcsos.x;
sos.g = qcsos.g;
sos.f = 0;
sos.p = [qcsos.p; dvar];

% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = struct('lin',Nl,'sos',Ns,'dsos',Nds,'sdsos',Nsds);
sosopt.Kc = struct('lin',Ml,'sos',Ms,'dsos',Mds,'sdsos',Msds);
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
