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
if isfield(opts.Kx,'l'), Nl = opts.Kx.l; else, Nl = 0; end
if isfield(opts.Kx,'s'), Ns = opts.Kx.s; else, Ns = 0; end
if isfield(opts.Kc,'l'), Ml = opts.Kc.l; else, Ml = 0; end
if isfield(opts.Kc,'s'), Ms = opts.Kc.s; else, Ms = 0; end

assert(n == (Nl + Ns), 'Dimension of Kx must be equal to number of variables (%d).', n);
assert(m == (Ml + Ms), 'Dimension of Kc must be equal to number of constraints (%d).', m)

% build SOS problem
sos.x = qcsos.x;
sos.g = qcsos.g;
sos.f = 0;
sos.p = [qcsos.p; dvar];

% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = struct('l',Nl,'s',Ns);
sosopt.Kc = struct('l',Ml,'s',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% store basis
% see SOSOPTCOMMON#GET_MONOMIALS_IN for details
obj.monom_xl = obj.sossolver.monom_xl;
obj.monom_xs = obj.sossolver.monom_xs;
obj.monom_p  = obj.sossolver.monom_p;
obj.monom_f  = obj.sossolver.monom_f;
obj.monom_gl = obj.sossolver.monom_gl;
obj.monom_gs = obj.sossolver.monom_gs;

end
