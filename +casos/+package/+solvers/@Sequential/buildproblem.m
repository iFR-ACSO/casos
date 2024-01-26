function buildproblem(obj,nlsos)
% Build SOS problem for bisection.

opts = obj.opts;


% problem size
n = length(nlsos.x);
m = length(nlsos.g);

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
sos.x = nlsos.x;
sos.g = nlsos.g;
sos.f = nlsos.f;
sos.p = nlsos.p;

% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = struct('l',Nl,'s',Ns);
sosopt.Kc = struct('l',Ml,'s',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% store basis
obj.monom_xl = basis(nlsos.x,~Is);
obj.monom_xs = basis(nlsos.x, Is);
obj.monom_p  = basis(nlsos.p);
obj.monom_f  = basis(nlsos.f);
obj.monom_gl = basis(nlsos.g,~Js);
obj.monom_gs = basis(nlsos.g, Js);

end
