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

% select sum-of-squares variables and constraints
Is = [false(Nl,1); true(Ns,1)];
Js = [false(Ml,1); true(Ms,1)];

% sparsity patterns
base_x = sparsity(nlsos.x);

% Jacobians
Df = jacobian(nlsos.f,nlsos.x);
Dg = jacobian(nlsos.g,nlsos.x);
% Hessian
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
sos_qp.derivatives.Df = Df;
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


%% constrained violation check via signed distance
% sparsity patterns
base_g = sparsity(nlsos.g);

r  = casos.PS.sym('r',length(nlsos.g));
s0 = casos.PD(base_g,ones(base_g.nnz,1));

sos_conVio = [];
% decision variables
sos_conVio.x = r;
% cost function
sos_conVio.f = sum(r);
% constraints
sos_conVio.g = nlsos.g + r.*s0;
% parameters to subproblem
% * parameters to nonlinear problem, current (nonlinear) solution
sos_conVio.p = [nlsos.p; nlsos.x];

% SOS options
sosopt               = opts.sossol_options;
sosopt.Kx.lin        = length(r);
sosopt.Kx.sos        = 0;
sosopt.Kc.sos        = length(sos_conVio.g);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.solver_conVio = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_conVio,sosopt);

%% Second-order correction (soc)

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
% augmented linear constraints (see Nocedal p.541)
dk        = x_star_par - nlsos.x ;
% actually xk cancelles out; just for readability 
sos_soc.g = dot(Dg,d) + conFun(nlsos.x + dk,nlsos.p) - dot(Dg, dk);
% parameters to subproblem
% * parameters to nonlinear problem, current (nonlinear) solution, Hessian,
% solution QP
sos_soc.p = [nlsos.p; nlsos.x; B(:);x_star_par];

sos_soc.derivatives.Hf = Hf;
sos_soc.derivatives.Df = Df;
sos_soc.derivatives.Dg = Dg;


% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = opts.Kx;
sosopt.Kc = opts.Kc;
sosopt.error_on_fail = false;

% initialize adapted convex SOS solver for soc
obj.solver_soc = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_soc,sosopt);

%% Damped BFGS (see Nocedal p.536/537)
lam_gs    = casos.PS.sym('lam_gs', obj.sparsity_gs);
x_k1      = casos.PS.sym('x_k1',base_x);

Langrangian = nlsos.f + dot(lam_gs,nlsos.g);

% search direction
obj.eval_s = casos.Function('f',{poly2basis(nlsos.x),poly2basis(x_k1),poly2basis(nlsos.p)}, { poly2basis( x_k1 - nlsos.x)});
% delta gradient of Langrangian 
obj.eval_y = casos.Function('f',{poly2basis(nlsos.x),poly2basis(x_k1),poly2basis(nlsos.p),poly2basis(lam_gs)}, {(op2basis(jacobian(Langrangian,x_k1)) - op2basis(jacobian(Langrangian,nlsos.x)) )' });

% using MX instead of SX is faster to build; no significant execution speed
% during online execution (but just a few examples)
theta = casadi.MX.sym('theta');
s     = casadi.MX.sym('s',[size(B,1),1]);
y     = casadi.MX.sym('y',[size(B,1),1]);
r     = casadi.MX.sym('r',[size(B,1),1]);
B     = casadi.MX.sym('B',[size(B,1),size(B,1)]);

% Powell-damping
obj.eval_r = casos.Function('f',{B,theta,y,s}, {theta*y+(1-theta)*B*s});

% B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s)
obj.damped_BFGS =  casadi.Function('f',{B,r,s}, {B - (B*(s*s')*B)/(s'*B*s) + (r*r')/(s'*r)});


%% feasibility restoration phase


sos_feas = [];
% convex decision variables
sos_feas.x = casos.PS.sym('x_feas',base_x);
% search direction
d = sos_feas.x - nlsos.x;
% quadratic cost function
sos_feas.f = (1/2)*dot2(Hf,d,d) + dot(d,d);
% linear constraints
sos_feas.g = nlsos.g + dot(Dg,d);
% parameters to subproblem
% * parameters to nonlinear problem, current (nonlinear) solution, Hessian
sos_feas.p = [nlsos.p; nlsos.x; B(:)];

% SOS options
sosopt               = opts.sossol_options;
sosopt.Kx            = opts.Kx;
sosopt.Kc            = opts.Kc;
sosopt.error_on_fail = false;

sos_feas.derivatives.Hf = Hf;
sos_feas.derivatives.Df = Df;
sos_feas.derivatives.Dg = Dg;

% initialize convex SOS solver
obj.solver_feasRes = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_qp,sosopt);


%% evaluation function
% nonlinear cost function 
obj.eval_cost      = casos.Function('f',{poly2basis(nlsos.x),poly2basis(nlsos.p)}, { full(nlsos.f) });

obj.eval_gradCost  = casos.Function('f',{poly2basis(nlsos.x),poly2basis(nlsos.p)}, { op2basis( jacobian(nlsos.f,nlsos.x) ) });

% gradient langrangian
obj.eval_gradLang =  casos.Function('f',{poly2basis(nlsos.x),poly2basis(nlsos.p),poly2basis(lam_gs)}, {Fnorm2( jacobian(Langrangian,nlsos.x) )' });


%% get parameter that are needed for initilization
obj.init_para.no_dual_var   = obj.sparsity_gs.nnz;
obj.init_para.size_B        = nnz(base_x);
obj.init_para.conVio.no_con = length(sos_conVio.g);


%% get params for display output (problem size etc)
obj.display_para.no_decVar         = obj.solver_convex.sparsity_x.nnz;
obj.display_para.no_sosCon         = length(nlsos.g);
obj.display_para.solver            = obj.opts.sossol;
obj.display_para.solver_build_time = toc(buildTime_in);


end