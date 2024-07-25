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

%% build SOS problem of underlying Q-SDP

% store user parameter
p0 = nlsos.p;

base_x = sparsity(nlsos.x);

% current iterate
xi_k    = casos.PS.sym('xi_k',base_x);

% search direction
d       = casos.PS.sym('d',base_x);

% search direction is decision variable of underlying Q-SDP
sos.x = d;

% Taylor Approximation constraints and evaluate at current solution
conFun      = casos.Function('f',{nlsos.x, p0},{nlsos.g});
derivConFun = casos.Function('f',{nlsos.x, d,p0},{ dot(jacobian(nlsos.g,nlsos.x),d) });

sos.g = conFun(xi_k,p0) + derivConFun(xi_k,d,p0); % sos.g = linearize(nlsos.g,nlsos.x,xi_k);

% parameterize cost in hessian
B_k    = casos.PS.sym('b',[length(poly2basis(xi_k)),length(poly2basis(xi_k))]);


sos.derivatives.Hf = casadi.SX(B_k);

% needed to initialize it later
obj.sizeHessian = size(B_k);

% needed because we have a parameter vector
B_k = B_k(:);


% quadratic cost approximation; refactor hessian vector to symmetric matrix
% f =          1/2 d^T B d + nabla f(x_k)^T*d
nabla_f_Fun = casos.Function('f',{nlsos.x,d,p0},{ dot(jacobian(nlsos.f,nlsos.x),d) });

sos.f =  1/2*poly2basis(d)'* reshape(B_k,[length(poly2basis(xi_k)),length(poly2basis(xi_k))]) * poly2basis(d) + nabla_f_Fun(xi_k,d,p0); %linearize(nlsos.f,sos.x,xi_k);

% extend parameter vector
sos.p = [p0; xi_k; B_k];

% initilize Filter object
obj.Filter = Filter([]);

% SOS options
sosopt               = opts.sossol_options;
sosopt.Kx            = struct('lin',Nl,'sos',Ns);
sosopt.Kc            = struct('lin',Ml,'sos',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver (subproblem)
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

%% Second-order correction

% correction term
correction = conFun(xi_k + d ,p0) - derivConFun(xi_k,d,p0);

% new decision variable is corrected search direction
dsoc = casos.PS.sym('dsoc',base_x);

% get adapted constraint
sosSOC.g = conFun(xi_k,p0) + derivConFun(xi_k,dsoc,p0) + correction; %linearize(nlsos.g,nlsos.x,xi_k) + correction;

sosSOC.p = [p0; xi_k; d; B_k];

sosSOC.x = dsoc;

sosSOC.f = 0;

% initialize SOS solver for SOC
obj.solver_soc = casos.package.solvers.sossolInternal('SOS',opts.sossol,sosSOC,sosopt);

%% setup damped BFGS
lam_gs    =  casos.PS.sym('lam_gs', obj.sparsity_gs);

s = casadi.MX.sym('s',[length(poly2basis(nlsos.x)),1]);
r = casadi.MX.sym('r',[length(poly2basis(nlsos.x)),1]);
B = casadi.MX.sym('B',[length(poly2basis(nlsos.x)),length(poly2basis(nlsos.x))]);

obj.BFGS_fun =  casadi.Function('f',{B,r,s}, {B + (r*r')/(s'*r) - (B*(s*s')*B)/(s'*B*s)});
 
% Langrangian
dLdx = jacobian(nlsos.f + dot(lam_gs,nlsos.g),nlsos.x)';

% gradient of Langrangian needed for BFGS and for convergence check
obj.nabla_xi_L      = casos.Function('f',{poly2basis(nlsos.x),poly2basis(lam_gs),p0}, { op2basis(dLdx) });
obj.nabla_xi_L_norm = casos.Function('f',{poly2basis(nlsos.x),poly2basis(lam_gs),p0}, { Fnorm2(dLdx) });


% cost function and gradient needed for filter linesearch
obj.f          = casos.Function('f',{poly2basis(nlsos.x),p0}, { nlsos.f });
obj.nabla_xi_f = casos.Function('f',{poly2basis(nlsos.x),p0}, { op2basis(jacobian(nlsos.f,nlsos.x)) });


%% setup projection for constraint violation check

% work around for polynomial interface
obj.xk1fun = casos.Function('f',{poly2basis(nlsos.x),p0}, {nlsos.x});

% identify nonlinear constraints 
I = zeros(length(nlsos.g),1);
for idx = 1:length(nlsos.g)
    I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
end

% Gram decision variable
s    = casos.PS.sym('q',sparsity(nlsos.g(I==1)));

if ~isempty(s)

    % parameterized projection for linesearch prediction
    nonLinCon = nlsos.g(I==1);
    conFunRed = casos.Function('f',{nlsos.x,p0}, {nonLinCon});
    
    % projection error
    e = s - conFunRed(nlsos.x,p0);
    
    % define Q-SOS problem parameterized nonlinear constraints
    %   min ||s-p||^2  s.t. s is SOS
    proj_sos = struct('x',s,'f',dot(e,e),'g',s,'p',[nlsos.x;p0]);
    
    opts               = [];
    opts.Kc            = struct('sos', length(s));
    opts.error_on_fail = 0;
    obj.projConPara    =  casos.sossol('S','mosek',proj_sos,opts);
else

    obj.projConPara    = [];

end


%% feasibility restoration phase
    
% get multiplier for SOS cone projection
s    = casos.PS.sym('q',grambasis(nlsos.g(I==1))); 

% Currently just debugging for feasibility restoration
x    = casos.Indeterminates('x',2,1);

obj.size_s = length(s);
obj.s0     = ones(obj.size_s,1)*(x'*x);
obj.size_x = length(sos.x);


% projection error
e = s - nlsos.g(I==1);

% weight for regularization (parameter)
zeta = casos.PS.sym('z');

% weight for regularization (parameter)
minCost = casos.PS.sym('m');

% current iterate where feas. restoration is called (parameter)
x0   = casos.PS.sym('xi0',base_x);


eReg = nlsos.x-x0;
regularization = dot(eReg,eReg);

cost = dot(e,e) + zeta/2*regularization;

sosFeas = struct('x',[nlsos.x; s],...      % augment decision variables
                 'f',cost , ...
                 'p',[x0;zeta;minCost]);

sosFeas.('g') = [nlsos.g(I~=1);s];

opts               = [];
opts.Kc            = struct('sos', length(sosFeas.g));
opts.Kx            = struct('lin', length(sosFeas.x));
opts.error_on_fail = 1;
opts.verbose       = 0;

% initialize solver
obj.solver_feas_res = casos.nlsossol('S','FeasRes',sosFeas,opts);

end
