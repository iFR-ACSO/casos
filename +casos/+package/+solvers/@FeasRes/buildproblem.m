function buildproblem(obj,nlsos)
% Build SOS problem for sequential SOS.
import casos.package.solvers.globalization.Filter
opts = obj.opts;


% problem size
n = length(nlsos.x);
m = length(nlsos.g);

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
sos.x = nlsos.x;
sos.g = nlsos.g;
sos.f = nlsos.f;
sos.p = nlsos.p;

% store user parameter
p0 = nlsos.p;


base_x = sparsity(sos.x);

xi_k    = casos.PS.sym('xi_k',base_x);
d       = casos.PS.sym('d',base_x);

% Taylor Approximation constraints and evaluate at current solution
sos.g = linearize(nlsos.g,nlsos.x,xi_k);

% parameterize cost in hessian
B_k    = casos.PS.sym('b',[length(poly2basis(xi_k)),length(poly2basis(xi_k))]);

obj.sizeHessian = size(B_k);
sos.derivatives.Hf = casadi.SX(B_k);

% needed because we have a parameter vector
B_k = B_k(:);

% gf_fun = casos.Function('f',{poly2basis(nlsos.x)},{op2basis(jacobian(nlsos.f,nlsos.x))});
% sos.derivatives.Gf = casadi.SX(gf_fun(sos.x));


% gg_fun = casos.Function('f',{nlsos.x},{op2basis(jacobian(nlsos.g,nlsos.x))});
% sos.derivatives.Gg = casadi.SX(gg_fun(sos.x));


% quadratic cost approximation; refactor hessian vector to symmetric matrix
sos.f =  1/2*poly2basis(sos.x)'* reshape(B_k,[length(poly2basis(xi_k)),length(poly2basis(xi_k))]) * poly2basis(sos.x) + linearize(nlsos.f,nlsos.x,xi_k);

% extend parameter vector
sos.p = [p0; xi_k; B_k];

% initilize Filter object
obj.Filter = Filter([]);

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
obj.sparsity_p  = sparsity(nlsos.p);
obj.sparsity_f  = sparsity(nlsos.f);
obj.sparsity_g  = obj.sossolver.sparsity_g;
obj.sparsity_gl = obj.sossolver.sparsity_gl;
obj.sparsity_gs = obj.sossolver.sparsity_gs;

%% setup damped BFGS
lam_gs    =  casos.PS.sym('lam_gs', obj.sparsity_gs);

s = casadi.MX.sym('s',[length(poly2basis(sos.x)),1]);
r = casadi.MX.sym('r',[length(poly2basis(sos.x)),1]);
B = casadi.MX.sym('B',[length(poly2basis(sos.x)),length(poly2basis(sos.x))]);

obj.BFGS_fun =  casadi.Function('f',{B,r,s}, {B + (r*r')/(s'*r) - (B*(s*s')*B)/(s'*B*s)});
 
% Langrangian
dLdx = jacobian(nlsos.f + dot(lam_gs,nlsos.g),sos.x)';

% gradient of Langrangian needed for BFGS and for convergence check
obj.nabla_xi_L      =  casos.Function('f',{poly2basis(sos.x),poly2basis(lam_gs),poly2basis(p0)}, { op2basis(dLdx) });
obj.nabla_xi_L_norm = casos.Function('f',{poly2basis(sos.x),poly2basis(lam_gs),poly2basis(p0)}, { Fnorm2(dLdx) });

% cost function and gradient needed for filter linesearch
obj.f          = casos.Function('f',{poly2basis(sos.x),poly2basis(p0)}, { nlsos.f });
obj.nabla_xi_f = casos.Function('f',{poly2basis(sos.x),poly2basis(p0)}, { op2basis(jacobian(nlsos.f,sos.x)) });


%% setup projection for constraint violation check

% % identify nonlinear constraints; work around for polynomial interface
obj.xk1fun = casos.Function('f',{poly2basis(sos.x),poly2basis(p0)}, {sos.x});

% error should be zero anyway
obj.projConPara    = [];


end
