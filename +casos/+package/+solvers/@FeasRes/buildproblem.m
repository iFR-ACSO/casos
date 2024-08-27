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

% distinguish between linear and sos decision variables
nonLinProb_linDec = nlsos.x(1:Nl);
nonLinProb_sosDec = nlsos.x(Nl+1:end);

base_x_lin = sparsity(nonLinProb_linDec);
base_x_sos = grambasis(nonLinProb_sosDec);

% current iterate (just for parameterization)
xi_k_lin    = casos.PS.sym('xi_kl',base_x_lin);
xi_k_sos    = casos.PS.sym('xi_ks',base_x_sos);

xi_k = [xi_k_lin;xi_k_sos];

% new convex solution
xi_new_lin = casos.PS.sym('xi_l',base_x_lin);
xi_new_sos = casos.PS.sym('xi_s',base_x_sos);

xi_new = [xi_new_lin;xi_new_sos];

% search direction
d = xi_new - xi_k;

d_new_lin = casos.PS.sym('d_l',base_x_lin);
d_new_sos = casos.PS.sym('d_s',base_x_sos);

d_star = [d_new_lin;d_new_sos];

% search direction is decision variable of underlying Q-SDP
sos.x = xi_new;

% Taylor Approximation constraints and evaluate at current solution
conFun       = casos.Function('f',{nlsos.x, p0},{nlsos.g});
dgdxi        = jacobian(nlsos.g,nlsos.x);
derivConFun  = casos.Function('f',{nlsos.x, xi_new,p0},{ dot(dgdxi,xi_new-nlsos.x) });
derivConFunC = casos.Function('f',{nlsos.x, d_star,p0},{ dot(dgdxi,d_star) });

% parameterize in d_star; we ony have to build one SDP solver for both
% convex subproblem and SOC
sos.g = conFun(xi_k+d_star,p0) + derivConFun(xi_k, xi_new,p0) -  derivConFunC(xi_k, d_star,p0);


% parameterize cost in hessian
size_bk = length(poly2basis(xi_k));

B_k    = casos.PS.sym('b',[size_bk,size_bk]);

sos.derivatives.Hf = casadi.SX(B_k);

% needed to initialize it later
obj.sizeHessian = size(B_k);


% quadratic cost approximation; refactor hessian vector to symmetric matrix
% f =          1/2 d^T B d + nabla f(x_k)^T*d
nabla_f_Fun = casos.Function('f',{nlsos.x, xi_new,p0},{ dot(jacobian(nlsos.f,nlsos.x),xi_new-nlsos.x) });

sos.f =  1/2*poly2basis(d)'* B_k * poly2basis(d) + nabla_f_Fun(xi_k, xi_new,p0); %linearize(nlsos.f,sos.x,xi_k);

% extend parameter vector
sos.p = [p0; xi_k; B_k(:);d_star];

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

%% setup damped BFGS
lam_gs    =  casos.PS.sym('lam_gs', obj.sparsity_gs);

s = casadi.MX.sym('s',[size_bk,1]);
r = casadi.MX.sym('r',[size_bk,1]);
B = casadi.MX.sym('B',[size_bk,size_bk]);

obj.BFGS_fun =  casadi.Function('f',{B,r,s}, {B + (r*r')/(s'*r) - (B*(s*s')*B)/(s'*B*s)});

%%  Pre-compute Langrangian derivative, cost derivative etc. (needed for filter and/or convergence check)
dLdx = jacobian(nlsos.f + dot(lam_gs,nlsos.g),nlsos.x)';

% gradient of Langrangian needed for BFGS and for convergence check
[coeff_nlsos_x,~] = poly2basis(nlsos.x);
coeff_lam_gs                 = poly2basis(lam_gs);
coeff_p0                     = poly2basis(p0);

obj.nabla_xi_L      = casos.Function('f',{coeff_nlsos_x,coeff_lam_gs,coeff_p0}, { op2basis(dLdx) });
obj.nabla_xi_L_norm = casos.Function('f',{coeff_nlsos_x,coeff_lam_gs,coeff_p0}, { Fnorm2(dLdx) });

% cost function and gradient needed for filter linesearch
obj.f          = casos.Function('f',{coeff_nlsos_x,coeff_p0}, { nlsos.f });
obj.nabla_xi_f = casos.Function('f',{coeff_nlsos_x,coeff_p0}, { op2basis(jacobian(nlsos.f,nlsos.x)) });

obj.nabla_xi_g = casos.Function('f',{coeff_nlsos_x,coeff_p0}, { op2basis(dgdxi) });


%% setup projection for constraint violation check

% work around for polynomial interface
obj.xk1fun = casos.Function('f1',{coeff_nlsos_x ,coeff_p0}, {nlsos.x});
obj.p0poly= casos.Function('f2',{coeff_p0}, {p0});

% identify nonlinear constraints
I = zeros(length(nlsos.g),1);
for idx = 1:length(nlsos.g)
    I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
end
    
    % decision variable (linear dec. var + sos constraint)
    s    = casos.PS.sym('q',grambasis(nlsos.g(I==1)));
    
    if ~isempty(s)
        
        % conFun = casos.Function('g1',{nlsos.x,p0},{nlsos.g(I==1)});
        
        % projection error
        e = s - nlsos.g(I==1);
        
        
        % define Q-SOS problem parameterized nonlinear constraints
        %   min ||s-p||^2  s.t. s is SOS
        proj_sos = struct('x',s,...
                          'g',s,...
                        'f',pnorm2(e), ...
                        'p',[nlsos.x;p0]);
        
        opts               = [];
        opts.Kx            = struct('lin', length(s));
        opts.error_on_fail = 1;
        obj.projConPara    =  casos.sossol('S','mosek',proj_sos,opts);
        
        
    else
        
        obj.projConPara    = [];
        
    end
     
    
end


