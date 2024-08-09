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

% distinguish between linear and sos dec. var
nonLinProb_linDec = nlsos.x(1:Nl);
nonLinProb_sosDec = nlsos.x(Nl+1:end);

base_x_lin = sparsity(nonLinProb_linDec);
base_x_sos = grambasis(nonLinProb_sosDec);

% current iterate
xi_k_lin    = casos.PS.sym('xi_kl',base_x_lin);
xi_k_sos    = casos.PS.sym('xi_ks',base_x_sos);

xi_k = [xi_k_lin;xi_k_sos];

% new convex solution
xi_new_lin = casos.PS.sym('xi_l',base_x_lin);
xi_new_sos = casos.PS.sym('xi_s',base_x_sos);

xi_new = [xi_new_lin;xi_new_sos];

d = xi_new - xi_k; 


d_new_lin = casos.PS.sym('d_l',base_x_lin);
d_new_sos = casos.PS.sym('d_s',base_x_sos);

d_star = [d_new_lin;d_new_sos];

% search direction is decision variable of underlying Q-SDP
sos.x = xi_new; 

% Taylor Approximation constraints and evaluate at current solution
conFun      = casos.Function('f',{nlsos.x, p0},{nlsos.g});
derivConFun = casos.Function('f',{nlsos.x, xi_new,p0},{ dot(jacobian(nlsos.g,nlsos.x),xi_new-nlsos.x) });
derivConFunC = casos.Function('f',{nlsos.x, d_star,p0},{ dot(jacobian(nlsos.g,nlsos.x),d_star) });

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
 
% Langrangian
dLdx = jacobian(nlsos.f + dot(lam_gs,nlsos.g),nlsos.x)';

% gradient of Langrangian needed for BFGS and for convergence check
[coeff_nlsos_x,sparse_nlsos] = poly2basis(nlsos.x);
coeff_lam_gs  = poly2basis(lam_gs);
coeff_p0      = poly2basis(p0);

obj.nabla_xi_L      = casos.Function('f',{coeff_nlsos_x,coeff_lam_gs,p0}, { op2basis(dLdx) });
obj.nabla_xi_L_norm = casos.Function('f',{coeff_nlsos_x,coeff_lam_gs,p0}, { Fnorm2(dLdx) });

% cost function and gradient needed for filter linesearch
obj.f          = casos.Function('f',{coeff_nlsos_x,p0}, { nlsos.f });
obj.nabla_xi_f = casos.Function('f',{coeff_nlsos_x,p0}, { op2basis(jacobian(nlsos.f,nlsos.x)) });

obj.nabla_xi_g = casos.Function('f',{coeff_nlsos_x,p0}, { op2basis(jacobian(nlsos.g,nlsos.x)) });


%% setup projection for constraint violation check

% work around for polynomial interface
obj.xk1fun = casos.Function('f1',{coeff_nlsos_x ,coeff_p0}, {nlsos.x});

% identify nonlinear constraints 
I = zeros(length(nlsos.g),1);
for idx = 1:length(nlsos.g)
    I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
end

if  strcmp(obj.opts.conViolCheck,'projection')
    
    
    % Gram decision variable
    s    = casos.PS.sym('q',grambasis(nlsos.g(I==1)));
    
    if ~isempty(s) 
        
        conFun = casos.Function('g1',{nlsos.x},{nlsos.g(I==1)});
        
        % projection error
        e = s - conFun(xi_k);
        
        
        % define Q-SOS problem parameterized nonlinear constraints
        %   min ||s-p||^2  s.t. s is SOS
        proj_sos = struct('x',s,'g',s,'f',dot(e,e),'p',xi_k);
        
        opts               = [];
        opts.Kx            = struct('lin', length(s));
        opts.Kc            = struct('sos', length(s));
        opts.error_on_fail = 0;
        obj.projConPara    =  casos.sossol('S','mosek',proj_sos,opts);
        

    else
    
        obj.projConPara    = [];
    
    end


else
    % constraint violation is checked via pseudo-proejction
    
    % did the user provide samples
    if isempty(opts.conVioSamp)

        % if no, we simply build a scaled (scaling set to 10) 
        % unit box for all indeterminates
        nIndet = length(sparse_nlsos.indeterminates);

        a = -10;
        b =  10;

        obj.opts.conVioSamp = a + (b-a)*rand(nIndet,1000);

    end

    nonLinCon = nlsos.g(I==1);
    
    % evaluate current solution (coefficients) at provided sampling points
    obj.pseudoProj = casos.Function('f',{coeff_nlsos_x,p0}, {nonLinCon});
    
    
    % for the final-check, we should still use the actual proejction to
    % check the result

    % Gram decision variable
    s    = casos.PS.sym('q',grambasis(nlsos.g(I==1)));
    
    if ~isempty(s) 
        
        conFun = casos.Function('g1',{nlsos.x},{nlsos.g(I==1)});
        
        % projection error
        e = s - conFun(xi_k);
        
        
        % define Q-SOS problem parameterized nonlinear constraints
        %   min ||s-p||^2  s.t. s is SOS
        proj_sos = struct('x',s,'g',s,'f',dot(e,e),'p',xi_k);
        
        opts               = [];
        opts.Kx            = struct('lin', length(s));
        opts.Kc            = struct('sos', length(s));
        opts.error_on_fail = 0;
        obj.projConPara    =  casos.sossol('S','mosek',proj_sos,opts);
        

    else
    
        obj.projConPara    = [];
    
    end



end


%% feasibility restoration phase

if obj.opts.feasRes_actv_flag
% get multiplier for SOS cone projection

spars_non_g = sparsity(nlsos.g(I==1));
s           = casos.PS.sym('q',spars_non_g ); 


obj.size_s = length(s);
obj.s0     = casos.PD(spars_non_g , ones(spars_non_g.nnz,1));
obj.size_x = length(sos.x);


% projection error
e = s - nlsos.g(I==1);

% weight for regularization (parameter)
zeta = casos.PS.sym('z');
vio0 = casos.PS.sym('V');


% current iterate where feas. restoration is called (parameter)

% current iterate
xi_0_lin    = casos.PS.sym('x0',base_x_lin);
xi_0_sos    = casos.PS.sym('x0',base_x_sos);

x0 = [xi_0_lin;xi_0_sos];

eReg           = nlsos.x - x0;

regularization = dot(eReg,eReg);

cost = 1/2*dot(e,e) + 1/2*regularization*zeta;


obj.conVio_0 = casos.Function('f',{nlsos.x,s,x0,zeta},{cost});

sosFeas = struct('x',[nlsos.x; s],...      % augment decision variables
                 'f',cost , ...
                 'p',[zeta;vio0;x0]);

sosFeas.('g') = [nlsos.g(I==1);s];

opts               = [];
opts.Kc            = struct('sos', length(sosFeas.g));
opts.Kx            = struct('lin', length(sosFeas.x));
opts.error_on_fail = 1;
opts.verbose       = 1;

% initialize solver
obj.solver_feas_res = casos.nlsossol('S','FeasRes',sosFeas,opts);

else
    % feasibility restoration is not active!
    obj.solver_feas_res = [];

end

end
