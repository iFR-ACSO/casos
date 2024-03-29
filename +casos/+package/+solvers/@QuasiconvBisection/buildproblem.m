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
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% store basis
% see SOSOPTCOMMON#GET_MONOMIALS_IN for details
% <<<<<<< HEAD
% obj.monom_xl = obj.sossolver.monom_xl;
% obj.monom_xs = obj.sossolver.monom_xs;
% obj.monom_p  = obj.sossolver.monom_p;
% obj.monom_f  = obj.sossolver.monom_f;
% obj.monom_gl = obj.sossolver.monom_gl;
% obj.monom_gs = obj.sossolver.monom_gs;
% =======
obj.monom_xl = monomials_in(obj.sossolver,2);
obj.monom_xs = monomials_in(obj.sossolver,4);
obj.monom_p  = basis(qcsos.p);
obj.monom_f  = basis(qcsos.f);
obj.monom_gl = monomials_in(obj.sossolver,5);
obj.monom_gs = monomials_in(obj.sossolver,7);
% >>>>>>> e46de63974e6a5759e850b41d75652d9cd4f52a5

end
