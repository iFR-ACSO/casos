function buildproblem(obj,nlsos)
% Build SOS problem for bisection.

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

% build SOS problem
sos.x = nlsos.x;
sos.g = nlsos.g;
sos.f = nlsos.f;
sos.p = nlsos.p;

% store user parameter
p0 = sos.p;
   
% if is_symbolic(sos.f)     
%     obj.signCost = 1;
% else
%     obj.signCost = -1;
% end


% for each decision variable generate a parameter (which is the current solution)
% naming is not important since this is only for internal 
base_x = basis(sos.x);

xi_k    = casos.PS.sym('xi_k',base_x);
xi_plus = casos.PS.sym('xi_plus',base_x);



% Taylor Approximation constraints and evaluate at current solution
sos.g = linearize(nlsos.g,nlsos.x,xi_k);

% Taylor Approximation cost and evaluate at current solution

% sos.f = linearize(nlsos.f,nlsos.x,xi_k);

% parameterize cost in hessian
H    = casos.PS.sym('h',[length(poly2basis(xi_k)),length(poly2basis(xi_k))]);

obj.sizeHess = size(H);
% needed because we have a parameter vector
H = H(:);

% quadratic cost approximation; refactor hessian vector to symmetric matrix
sos.f =  1/2*poly2basis(sos.x)'* reshape(H,[length(poly2basis(xi_k)),length(poly2basis(xi_k))]) *poly2basis(sos.x) + jacobian(sos.f,sos.x)*poly2basis(sos.x)

% extend parameter vector
sos.p = [p0; xi_k;H];

% casos function for nonlinear cost and nonlinear constraints
obj.cost_fun      = casos.Function('f',{sos.x}, {nlsos.f});
obj.constraintFun = casos.Function('f',{sos.x,p0}, {nlsos.g});

% SOS options
sosopt    = opts.sossol_options;
sosopt.Kx = struct('lin',Nl,'sos',Ns);
sosopt.Kc = struct('lin',Ml,'sos',Ms);

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% see SOSOPTCOMMON#GET_MONOMIALS_IN for details
obj.monom_xl = monomials_in(obj.sossolver,2);
obj.monom_xs = monomials_in(obj.sossolver,4);
obj.monom_p  = basis(nlsos.p);
obj.monom_f  = basis(nlsos.f);
obj.monom_gl = monomials_in(obj.sossolver,5);
obj.monom_gs = monomials_in(obj.sossolver,7);


%% setup Merit-function
lam_gs    =  casos.PS.sym('lam_gs',obj.monom_gs);
dual_k    =  casos.PS.sym('lam_gs',obj.monom_gs);
dual_plus =  casos.PS.sym('lam_gs',obj.monom_gs);

% Langrange function
L = casos.Function('f',{sos.x,lam_gs,p0}, {nlsos.f - dot(lam_gs,nlsos.g)});

obj.Merit = L;

obj.dLdx = casos.Function('f',{sos.x,lam_gs,p0}, { jacobian(nlsos.f + dot(lam_gs,nlsos.g),sos.x)'  });

obj.langrangeLinear = casos.Function('f',{sos.x,xi_k,p0,lam_gs,H}, {jacobian(sos.f + dot(lam_gs,sos.g),sos.x)'}); 

%% setup line search
d = casos.PS.sym('d');

Psi_d = L(xi_plus*d + (1-d)*xi_k , lam_gs,p0) - 1e-15*d^2; % see bisos implementation

%% setup projection
% identify nonlinear constraints

idx = 1:length(nlsos.g);
I = arrayfun(@(i) ~is_linear(nlsos.g,nlsos.x,idx==i), idx);

% Gram decision variable
s    = casos.PS.sym('q',basis(nlsos.g,I));
gsym = casos.PS.sym('gsym',basis(nlsos.g,I));

% store indices of nonlinear constraints
obj.idxNonlinCon = I;

% projection error
e = s - gsym;

% define Q-SOS problem parameterized nonlinear constraints
%   min ||s-p||^2  s.t. s is SOS
proj_sos = struct('x',s,'f',dot(e,e),'g',s,'p',gsym);

opts                = [];
opts.Kc             = struct('sos', length(s));
opts.error_on_fail  = 0;
obj.projCon         =  casos.sossol('S','mosek',proj_sos,opts);


% parameterized projection for linesearch prediction
nonLinCon     = nlsos.g(I);
obj.conFunRed = casos.Function('f',{sos.x,p0}, {nonLinCon});

% projection error
e = s - obj.conFunRed(sos.x,p0);

% define Q-SOS problem parameterized nonlinear constraints
%   min ||s-p||^2  s.t. s is SOS
proj_sos = struct('x',s,'f',dot(e,e),'g',s,'p',[sos.x;p0]);

opts               = [];
opts.Kc            = struct('sos', length(s));
opts.error_on_fail = 0;
obj.projConPara    =  casos.sossol('S','mosek',proj_sos,opts);


%% work around for efficient computations in sequential-call
obj.plusFun        = casos.Function('f',{sos.x,xi_k }, {sos.x + xi_k });

obj.vertcatFun     = casos.Function('f',{p0,xi_k }, {[p0;xi_k]});

% obj.norm2FunOptVar = casos.Function('f',{xi_k,xi_plus,dual_k, dual_plus}, ...
%                      { dot(xi_k,xi_plus),dot(dual_k,dual_plus) });

obj.norm2FunOptVar = casos.Function('f',{xi_k,xi_plus,dual_k, dual_plus}, ...
                     { pnorm2(xi_k-xi_plus),pnorm2(dual_k-dual_plus) });

glin_sol         = casos.PS.sym('glin',basis(sos.g));

obj.norm2FunVio  = casos.Function('f',{ xi_plus, p0, glin_sol }, ...
                                      { pnorm2( obj.constraintFun(xi_plus,p0) - glin_sol) } );


dopt = casos.PS.sym('d');
obj.updateLineSearch  = casos.Function('f',{xi_k, xi_plus, dual_k, dual_plus, dopt }, ...
                        {dopt*xi_plus    + (1-dopt)*xi_k, dopt*dual_plus  + (1-dopt)*dual_k});




obj.deltaOptVar =  casos.Function('f',{xi_k, xi_plus, dual_k, dual_plus}, ...
                   {xi_plus-xi_k, dual_plus-dual_k});


end
