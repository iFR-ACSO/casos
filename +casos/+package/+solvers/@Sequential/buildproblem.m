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

%% Second-order correction


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% We only have one SDP solver left which is parameterized in the search
% direction!!!!!
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% xisoc_new_lin = casos.PS.sym('xi_ls',base_x_lin);
% xisoc_new_sos = casos.PS.sym('xi_ss',base_x_sos);
%
% xisoc_new = [xisoc_new_lin;xisoc_new_sos];
%
% dsoc = xisoc_new - xi_k;
%
% % correction term
% correction = conFun(xi_k + xi_new ,p0) - derivConFun(xi_k, xi_new,p0);
%
% obj.correctFun = casos.Function('fc',{poly2basis(xi_k), poly2basis(xi_new),poly2basis(p0)},{correction});
%
% % % get adapted constraint
% sosSOC.g = derivConFun(xi_k,xisoc_new,p0) + correction;
%
% sosSOC.p = [p0; xi_k; xi_new; B_k(:)];
% %
% sosSOC.x = xisoc_new;

% sosSOC.f =  1/2*poly2basis(dsoc)'* B_k * poly2basis(dsoc) + nabla_f_Fun(xi_k, xisoc_new,p0); %linearize(nlsos.f,sos.x,xi_k);

% sosSOC.derivatives.Hf = casadi.SX(reshape(B_k,[size_bk,size_bk]));

% % initialize SOS solver for SOC
% obj.solver_soc = casos.package.solvers.sossolInternal('SOS',opts.sossol,sosSOC,sosopt);

%% setup damped BFGS
lam_gs    =  casos.PS.sym('lam_gs', obj.sparsity_gs);

s = casadi.MX.sym('s',[size_bk,1]);
r = casadi.MX.sym('r',[size_bk,1]);
B = casadi.MX.sym('B',[size_bk,size_bk]);

obj.BFGS_fun =  casadi.Function('f',{B,r,s}, {B + (r*r')/(s'*r) - (B*(s*s')*B)/(s'*B*s)});

%%  Pre-compute Langrangian derivative, cost derivative etc. (needed for filter and/or convergence check)

Langrangian = nlsos.f + dot(lam_gs,nlsos.g);

% first order optimality condition
dLdx = jacobian(Langrangian,nlsos.x)';

% gradient of Langrangian needed for BFGS and for convergence check
[coeff_nlsos_x,sparse_nlsos] = poly2basis(nlsos.x);
coeff_lam_gs                 = poly2basis(lam_gs);
coeff_p0                     = poly2basis(p0);

obj.nabla_xi_L      = casos.Function('f',{coeff_nlsos_x,coeff_lam_gs,coeff_p0}, { op2basis(dLdx) });
obj.nabla_xi_L_norm = casos.Function('f',{coeff_nlsos_x,coeff_lam_gs,coeff_p0}, { norm(op2basis(dLdx),inf)   }); %Fnorm2(dLdx)

% cost function and gradient needed for filter linesearch
obj.f          = casos.Function('f',{coeff_nlsos_x,coeff_p0}, { nlsos.f });
obj.nabla_xi_f = casos.Function('f',{coeff_nlsos_x,coeff_p0}, { op2basis(jacobian(nlsos.f,nlsos.x)) });
obj.nabla_xi_g = casos.Function('f',{coeff_nlsos_x,coeff_p0}, { op2basis(dgdxi) });


%% setup projection for constraint violation check

% work around for polynomial interface
obj.xk1fun = casos.Function('f1',{coeff_nlsos_x ,coeff_p0}, {nlsos.x});
obj.p0poly = casos.Function('f2',{coeff_p0}, {p0});

% identify nonlinear constraints
I = zeros(length(nlsos.g),1);
for idx = 1:length(nlsos.g)
    I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
end

% check if we use projection or pseudo-projection
if  strcmp(obj.opts.conViolCheck,'projection')
    
    
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
                        'f',dot(e,e), ...
                        'p',[nlsos.x;p0]);
        
        opts_proj               = [];
        opts_proj.Kx            = struct('lin', length(s));
        opts_proj.Kc            = struct('sos', length(s));
        opts_proj.error_on_fail = 1;
        obj.projConPara    =  casos.sossol('S','scs',proj_sos,opts_proj);
        
        
    else
        
        obj.projConPara    = [];
        
    end
    
    
else
    % constraint violation is checked via pseudo-projection i.e. sampling
    
    % did the user provide samples
    if isempty(opts.conVioSamp)
        
        % if no, we simply build a scaled (scaling set to 10)
        % unit box for all indeterminates
        nIndet = length(sparse_nlsos.indeterminates);
        
        a = -10;
        b =  10;
        
        obj.opts.conVioSamp = a + (b-a)*rand(nIndet,1000);
        
    end

    % evaluate current solution (coefficients) at provided sampling points
    obj.pseudoProj = casos.Function('f',{coeff_nlsos_x,coeff_p0}, {nlsos.g(I==1)});
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % for the final-check, we should still use the actual projection to
    % check the result!
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % Gram decision variable
    s    = casos.PS.sym('q',grambasis(nlsos.g(I==1)));
    
    if ~isempty(s)
        
        % conFun = casos.Function('g1',{nlsos.x,p0},{nlsos.g(I==1)});
        
        % projection error
        e = s - nlsos.g(I==1);
        
        
        % define Q-SOS problem parameterized nonlinear constraints
        %   min ||s-p||^2  s.t. s is SOS
        proj_sos = struct('x',s,'g',s,'f',dot(e,e),'p',[nlsos.x;p0]);
        
        opts_proj               = [];
        opts_proj.Kx            = struct('lin', length(s));
        opts_proj.Kc            = struct('sos', length(s));
        opts_proj.error_on_fail = 1;
        obj.projConPara    =  casos.sossol('S','scs',proj_sos,opts_proj);
        
        
    else
        
        obj.projConPara    = [];
        
    end
   
end


%% feasibility restoration phase

% Here, we only setup the adjusted underlying QP. In eval_on_basis we setup
% the actual restoration phase

if ~obj.opts.feasRes_actv_flag
    % feasibility restoration is not active!
    obj.solver_feas_res = [];   
else
        
% we have an augmented decision variable vector

% get projection multiplier of ALL constraints
s_k_feas      = casos.PS.sym('q',grambasis(nlsos.g));
s_new_feas    = casos.PS.sym('q',grambasis(nlsos.g));

% distinguish between linear and sos decision variables
nonLinProb_linDe_feasc = nlsos.x(1:Nl);
nonLinProb_sosDec_feas = nlsos.x(Nl+1:end);

base_x_lin = sparsity(nonLinProb_linDe_feasc);
base_x_sos = grambasis(nonLinProb_sosDec_feas);

% current iterate (just for parameterization)
xi_k_lin_feas    = casos.PS.sym('xi_kl',base_x_lin);
xi_k_sos_feas    = casos.PS.sym('xi_ks',base_x_sos);

xi_k_feas = [xi_k_lin_feas;xi_k_sos_feas;s_k_feas];

lenOrigDecVar = length([xi_k_lin_feas;xi_k_sos_feas]);

% parameter
xi_0_lin_feas    = casos.PS.sym('xi_kl',base_x_lin);
xi_0_sos_feas    = casos.PS.sym('xi_ks',base_x_sos);

xi_0_feas = [xi_0_lin_feas;xi_0_sos_feas];

% new convex solution
xi_new_lin_feas = casos.PS.sym('xi_l',base_x_lin);
xi_new_sos_feas = casos.PS.sym('xi_s',base_x_sos);

xi_new_feas = [xi_new_lin_feas;xi_new_sos_feas;s_new_feas];

% search direction
d_feas = xi_new_feas - xi_k_feas;

d_new_lin_feas = casos.PS.sym('d_l',base_x_lin);
d_new_sos_feas = casos.PS.sym('d_s',base_x_sos);

d_star_feas = [d_new_lin_feas;d_new_sos_feas];

% search direction is decision variable of underlying Q-SDP
sos_feas.x = xi_new_feas;


% adjusted cost function
zeta = casos.PS.sym('zeta');
conFun = casos.Function('g1',{nlsos.x,p0},{nlsos.g});

e1 = s_k_feas - conFun(xi_k_feas(1:lenOrigDecVar),p0);
e2 = xi_k_feas(1:lenOrigDecVar) - xi_0_feas;

feasResCost = dot(e1,e1) + zeta*dot(e2,e2);


% parameterize cost in hessian
size_rk = length(poly2basis(xi_k_feas));

R_k    = casos.PS.sym('b',[size_rk,size_rk]);

sos_feas.derivatives.Hf = casadi.SX(R_k);

% (nabla_xi_k_feas f_feas)^T * (d_feas)
% we have several constant terms: zeta and xi_0_feas i.e. iterate where
% restoration is invoked
nabla_f_Fun = casos.Function('f',{xi_k_feas, xi_new_feas, p0, xi_0_feas, zeta},{ dot(jacobian(feasResCost,xi_k_feas),xi_new_feas-xi_k_feas) });

sos_feas.f =  1/2*poly2basis(d_feas)'* R_k * poly2basis(d_feas) + nabla_f_Fun(xi_k_feas, xi_new_feas,p0,xi_0_feas,zeta); %linearize(nlsos.f,sos.x,xi_k);

 

sos_feas.p = [p0; xi_k_feas;xi_0_feas;zeta;R_k(:);d_star_feas];
sos_feas.g = [];

% SOS options
sosopt_feas               = opts.sossol_options;
sosopt_feas.Kx            = struct('lin',length(xi_k_lin_feas),'sos',length(xi_k_sos_feas) + length(s_k_feas) );
sosopt_feas.Kc            = struct('lin',0,'sos',length(sos_feas.g));
sosopt_feas.error_on_fail = false;

% initialize convex SOS solver (subproblem)
obj.sossolver_feas = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_feas,sosopt_feas);




end

end
