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
sosopt.Kx = struct('l',Nl,'s',Ns);
sosopt.Kc = struct('l',Ml,'s',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% store basis
obj.monom_xl = basis(qcsos.x,~Is);
obj.monom_xs = basis(qcsos.x, Is);
obj.monom_p  = basis(qcsos.p);
obj.monom_f  = basis(qcsos.f);
obj.monom_gl = basis(qcsos.g,~Js);
obj.monom_gs = basis(qcsos.g, Js);

end
