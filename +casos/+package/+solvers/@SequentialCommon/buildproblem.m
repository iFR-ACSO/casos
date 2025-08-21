function buildproblem(obj,nlsos)
% Build convex subprobems for sequential algorithms.
buildTime_in = tic;

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


% sparsity patterns
base_x = sparsity(nlsos.x);

% Jacobians
Df = jacobian(nlsos.f,nlsos.x);
Dg = jacobian(nlsos.g,nlsos.x);

% Hessian /BFGS matrix
B = casadi.SX.sym('B',nnz(base_x),nnz(base_x));

Hf = casos.PSOperator(B,base_x,base_x);

%% build quadratic SOS problem
sos_qp = [];
% convex decision variables
sos_qp.x = casos.PS.sym('x_convex',base_x);
% search direction
d = sos_qp.x - nlsos.x;
% quadratic cost function
sos_qp.f = (1/2)*dot2(Hf,d,d) + dot(Df,d);

% linear constraints
sos_qp.g = nlsos.g + dot(Dg,d);

% parameters to subproblem
% * parameters to nonlinear problem, current (nonlinear) solution, Hessian
sos_qp.p = [nlsos.p; nlsos.x; B(:)];

% SOS options
sosopt               = opts.sossol_options;
sosopt.Kx            = opts.Kx;
sosopt.Kc            = opts.Kc;
sosopt.error_on_fail = false;

sos_qp.derivatives.Hf = Hf;
sos_qp.derivatives.Df = Df + dot(casos.PSOperator.from_primal(d),Hf); % correct for shift in variables
sos_qp.derivatives.Dg = Dg;

% initialize convex SOS solver
obj.solver_convex = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_qp,sosopt);

% store basis
obj.sparsity_x  = obj.solver_convex.sparsity_x;
obj.sparsity_xl = obj.solver_convex.sparsity_xl;
obj.sparsity_xs = obj.solver_convex.sparsity_xs;
obj.sparsity_p  = sparsity(nlsos.p);
obj.sparsity_f  = sparsity(nlsos.f);
obj.sparsity_g  = obj.solver_convex.sparsity_g;
obj.sparsity_gl = obj.solver_convex.sparsity_gl;
obj.sparsity_gs = obj.solver_convex.sparsity_gs;

% store basis in a public struct
obj.sparsity_pat_para.sparsity_x  = obj.solver_convex.sparsity_x;
obj.sparsity_pat_para.sparsity_xl = obj.solver_convex.sparsity_xl;
obj.sparsity_pat_para.sparsity_xs = obj.solver_convex.sparsity_xs;
obj.sparsity_pat_para.sparsity_p  = sparsity(nlsos.p);
obj.sparsity_pat_para.sparsity_f  = sparsity(nlsos.f);
obj.sparsity_pat_para.sparsity_g  = obj.solver_convex.sparsity_g;
obj.sparsity_pat_para.sparsity_gl = obj.solver_convex.sparsity_gl;
obj.sparsity_pat_para.sparsity_gs = obj.solver_convex.sparsity_gs;


if strcmp(obj.opts.conVioCheck,'signed-distance')
    %% constrained violation check via signed distance
    % sparsity patterns of nonlinear constraints

    % check only the nonlinear constraints
    I = true(length(nlsos.g),1);

    % for idx = 1:length(nlsos.g)
    %     I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
    % end
    % get gram half-basis for nonlinear constraints
    [~,~,z] = grambasis(nlsos.g,I);

    % build unit vectors
    base_s0 = gramunit(z);

    r  = casos.PS.sym('r',sum(I));
    s0 = casos.PD(base_s0);

    sos_conVio = [];
    % decision variables
    sos_conVio.x = r;
    % cost function
    sos_conVio.f = sum(r);
    % constraints
    sos_conVio.g = nlsos.g(I) + r.*s0;
    % parameters to subproblem
    % * parameters to nonlinear problem, current (nonlinear) solution
    sos_conVio.p = [nlsos.p; nlsos.x];

    % SOS options
    sosopt               = opts.sossol_options;
    sosopt.Kx.lin        = length(r);
    sosopt.Kx.sos        = 0;
    sosopt.Kc.sos        = length(sos_conVio.g);
    sosopt.error_on_fail = true;

    % initialize convex SOS solver
    obj.solver_conVio = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_conVio,sosopt);
else

    if ~isempty(obj.opts.userSample)
        x_sample = obj.opts.userSample;
    else
        x_sample = [];
    end

    X  = casadi.SX.sym('x',nlsos.g.nvars,1);

    G = subs(nlsos.g, nlsos.g.indeterminates, X);

    obj.conFun = casadi.Function('g', {poly2basis(nlsos.x) X}, {casadi.SX(G)});

    obj.x_sample = x_sample;

end
%% Second-order correction (soc)
if obj.opts.Soc_is_enabled
    % parameter of current qp solution
    x_star_par   = casos.PS.sym('x_star',base_x);
    % constraint function to be evaluated at xk+d
    conFun       = casos.Function('f',{nlsos.x, nlsos.p},{nlsos.g});

    sos_soc = [];
    % convex decision variables
    sos_soc.x = casos.PS.sym('x_soc',base_x);
    % search direction
    d = sos_soc.x - nlsos.x;
    % quadratic cost function
    sos_soc.f = (1/2)*dot2(Hf,d,d) + dot(Df,d);
    % augmented linear constraints (see Nocedal p.544)
    dk        = x_star_par - nlsos.x ;
    % actually xk cancelles out; just for readability
    sos_soc.g = dot(Dg,d) + conFun(nlsos.x + dk,nlsos.p) - dot(Dg, dk);
    % parameters to subproblem
    % * parameters to nonlinear problem, current (nonlinear) solution, Hessian,
    % solution QP
    sos_soc.p = [nlsos.p; nlsos.x; B(:);x_star_par];

    sos_soc.derivatives.Hf = Hf;
    sos_soc.derivatives.Df = Df + dot(casos.PSOperator.from_primal(d),Hf);
    sos_soc.derivatives.Dg = Dg;

    % SOS options
    sosopt    = opts.sossol_options;
    sosopt.Kx = opts.Kx;
    sosopt.Kc = opts.Kc;
    sosopt.error_on_fail = false;

    % initialize adapted convex SOS solver for soc
    obj.solver_soc = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_soc,sosopt);
end

%% Damped BFGS (see Nocedal p.536/537)
lam_gs    = casos.PS.sym('lam_gs', obj.sparsity_g);
x_k1      = casos.PS.sym('x_k1',base_x);

Langrangian = nlsos.f + dot(lam_gs,nlsos.g);

basis_nlsos_x = poly2basis(nlsos.x);
basis_xk1     = poly2basis(x_k1);
basis_nlsos_p = poly2basis(nlsos.p);
basis_lam_gs  = poly2basis(lam_gs);

% Langrangian casos function for evaluation
obj.L = casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, {Langrangian });

dLdx = op2basis(jacobian(Langrangian,nlsos.x))';

% Hessian of Langrangian; in case we use BFGS, we do not have to explictly
% compute the Hessian
if ~strcmp(obj.opts.hessian_approx,'BFGS')
    L_xx         = op2basis(hessian(Langrangian,nlsos.x));
    obj.hess_fun = casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { L_xx });
end
% search direction
obj.eval_s = casos.Function('f',{basis_nlsos_x,basis_xk1,basis_nlsos_p}, { poly2basis( x_k1 - nlsos.x)});

% delta gradient of Langrangian
dLdx = casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, {dLdx });
obj.dLdx = dLdx;

% y =  \nabla_x L (x_k1,lambda_k1) - \nabla_x L (x_k,lambda_k1)
obj.eval_y = casos.Function('f',{basis_nlsos_x,basis_xk1,basis_nlsos_p,basis_lam_gs}, {dLdx(basis_xk1,basis_nlsos_p,basis_lam_gs) - dLdx(basis_nlsos_x,basis_nlsos_p,basis_lam_gs) });

% using MX instead of SX is faster to build; no significant difference in
% execution time during online execution (but just a few examples)
theta = casadi.MX.sym('theta');
s     = casadi.MX.sym('s',[size(B,1),1]);
y     = casadi.MX.sym('y',[size(B,1),1]);
r     = casadi.MX.sym('r',[size(B,1),1]);
B     = casadi.MX.sym('B',[size(B,1),size(B,1)]);

% Powell-damping
obj.eval_r = casos.Function('f',{B,theta,y,s}, {theta*y+(1-theta)*B*s});

% B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s)
obj.damped_BFGS =  casadi.Function('f',{B,r,s}, {B - (B*(s*s')*B)/(s'*B*s) + (r*r')/(s'*r)});

%% evaluation functions
% nonlinear cost function
obj.eval_cost      = casos.Function('f',{basis_nlsos_x,basis_nlsos_p}, { full(nlsos.f) });
% gradient cost
obj.eval_gradCost  = casos.Function('f',{basis_nlsos_x,basis_nlsos_p}, { op2basis( jacobian(nlsos.f,nlsos.x) ) });

% norm gradient langrangian
%obj.eval_gradLang =  casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { Fnorm2( jacobian(Langrangian,nlsos.x) )' });

obj.eval_gradLang =  casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { norm( op2basis(jacobian(Langrangian,nlsos.x)),inf ) });


obj.eval_gradLang2  =  casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { norm(dot(lam_gs,nlsos.g),inf) });

obj.eval_constraint_fun = casos.Function('f',{basis_nlsos_x,basis_nlsos_p}, {poly2basis(nlsos.g)});

%% get parameter that are needed for initilization
obj.init_para.no_dual_var   = obj.sparsity_g.nnz;
obj.init_para.size_B        = nnz(base_x);

if strcmp(obj.opts.conVioCheck,'signed-distance')
    obj.init_para.conVio.no_con = length(sos_conVio.g);
end

%% get params for display output (problem size etc)
obj.display_para.no_decVar         = obj.solver_convex.sparsity_x.nnz;
obj.display_para.no_sosCon         = length(nlsos.g);
obj.display_para.solver            = obj.opts.sossol;
obj.display_para.solver_build_time = toc(buildTime_in);

end
