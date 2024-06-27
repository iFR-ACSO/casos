function buildproblem(obj,nlsos)
% Build SOS problem for bisection.

opts = obj.opts;

feasRes = opts.FeasRes;
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
   
% for each decision variable generate a parameter (which is the current solution)
% naming is not important since this is only for internal 
base_x = basis(sos.x);

obj.base_x = base_x;
xi_k    = casos.PS.sym('xi_k',base_x);
xi_plus = casos.PS.sym('xi_plus',base_x);
searchdir = casos.PS.sym('xi_plus',base_x);
% Taylor Approximation constraints and evaluate at current solution
sos.g = linearize(nlsos.g,nlsos.x,xi_k);

% Taylor Approximation cost and evaluate at current solution
if strcmp(obj.opts.Sequential_Algorithm,'SLP')

    sos.f = linearize(nlsos.f,nlsos.x,xi_k);
    
    % extend parameter vector
    sos.p = [p0; xi_k];
     

    obj.sizeHess =  0;
else

% parameterize cost in hessian
H    = casos.PS.sym('h',[length(poly2basis(xi_k)),length(poly2basis(xi_k))]);

obj.sizeHess = size(H);

% needed because we have a parameter vector
H = H(:);

% quadratic cost approximation; refactor hessian vector to symmetric matrix
sos.f =  1/2*poly2basis(sos.x)'* reshape(H,[length(poly2basis(xi_k)),length(poly2basis(xi_k))]) * poly2basis(sos.x) + linearize(nlsos.f,nlsos.x,xi_k);

% extend parameter vector
sos.p = [p0; xi_k;H];

end

% casos function for nonlinear cost and nonlinear constraints
obj.cost_fun      = casos.Function('f',{sos.x},    {nlsos.f});
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


%% setup adapted solver for second order correction (SOC) i.e. solve the adapted convex subproblem
obj.samplingPoint_proj = opts.Sampling_Point;

obj.constraintFun = casos.Function('f',{sos.x,p0},{nlsos.g});

helpFun1 = casos.Function('f',{xi_k,p0},{nlsos.g});

% Jacbian at current iterate
[G,zi,zo] = jacobian(nlsos.g,xi_k);

% Constraints with second order correction
correction = helpFun1(xi_k + searchdir ,p0) - casos.PS(zo,G*poly2basis(searchdir,zi));
sos.g = linearize(nlsos.g,nlsos.x,xi_k) + correction;

sos.p = [p0; xi_k;searchdir;H];

% initialize SOS solver for SOC
obj.solverSOC = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);


%% setup Merit-function
lam_gs    =  casos.PS.sym('lam_gs',obj.monom_gs);
dual_k    =  casos.PS.sym('lam_gs',obj.monom_gs);
dual_plus =  casos.PS.sym('lam_gs',obj.monom_gs);


if ~strcmp(obj.opts.Sequential_Algorithm,'SLP')

     s = casadi.MX.sym('s',[length(poly2basis(sos.x)),1]);
     r = casadi.MX.sym('r',[length(poly2basis(sos.x)),1]);
     B = casadi.MX.sym('B',[length(poly2basis(sos.x)),length(poly2basis(sos.x))]);

     obj.BFGS_fun =  casadi.Function('f',{B,r,s}, {B + (r*r')/(s'*r) - (B*(s*s')*B)/(s'*B*s)});

end


% Langrange function
L = casos.Function('f',{sos.x,lam_gs,p0}, {nlsos.f + dot(lam_gs,nlsos.g)});

obj.Merit = L;

obj.dLdx = casos.Function('f',{sos.x,lam_gs,p0}, { jacobian(nlsos.f + dot(lam_gs,nlsos.g),sos.x)'  });

obj.nabla_f = casos.Function('f',{sos.x,p0}, { jacobian(nlsos.f,sos.x)' });

obj.cost = casos.Function('f',{sos.x,p0}, { nlsos.f });


%% setup projection for constraint violation check
% identify nonlinear constraints
idx = 1:length(nlsos.g);
I = arrayfun(@(i) ~is_linear(nlsos.g,nlsos.x,idx==i), idx);

% Gram decision variable
s    = casos.PS.sym('q',basis(nlsos.g,I));

% store indices of nonlinear constraints
obj.idxNonlinCon = I;

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


%% setup feasibility restoration
% projection error
if ~feasRes
e = s - obj.conFunRed(sos.x,p0);

% define Q-SOS problem parameterized nonlinear constraints
%   min ||s-p||^2  s.t. s is SOS
feasRes_sos = struct('x',[s;sos.x],'f',dot(e,e),'g',[s;sos.x],'p',p0);


optsFeas.Sequential_Algorithm = 'SQP';
optsFeas.Kx      = struct('lin', length(feasRes_sos.x));
optsFeas.Kc      = struct('sos', length(feasRes_sos.g));
optsFeas.verbose = 1;
optsFeas.sossol_options.sdpsol_options.error_on_fail = 0;
optsFeas.FeasRes = 1;

obj.FeasResPhase = casos.nlsossol('S1','sequential',feasRes_sos,optsFeas);


end

%% work around for efficient computations in sequential-call; polynomial operations "outsourced" to functions
obj.plusFun        = casos.Function('f',{sos.x,xi_k }, {sos.x + xi_k });

obj.vertcatFun     = casos.Function('f',{p0,xi_k }, {[p0;xi_k]});

% 
% obj.norm2FunOptVar = casos.Function('f',{xi_k,xi_plus,dual_k, dual_plus}, ...
%                      { pnorm2(xi_k-xi_plus), pnorm2(dual_k-dual_plus) });

obj.norm2FunOptVar = casos.Function('f',{xi_k,xi_plus,dual_k, dual_plus}, ...
                     { dot(xi_k-xi_plus,xi_k-xi_plus), dot(dual_k-dual_plus,dual_k-dual_plus) });


dopt = casos.PS.sym('d');
obj.updateLineSearch  = casos.Function('f',{xi_k, xi_plus, dual_k, dual_plus, dopt }, ...
                        {dopt*xi_plus    + (1-dopt)*xi_k, dopt*dual_plus  + (1-dopt)*dual_k});




obj.deltaOptVar =  casos.Function('f',{xi_k, xi_plus, dual_k, dual_plus}, ...
                                      {xi_plus-xi_k, dual_plus-dual_k});



obj.delta_search = casos.Function('f',{xi_plus,xi_k}, {poly2basis(xi_plus - xi_k,base_x)});


end
