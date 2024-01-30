function buildproblem(obj,nlsos)
% Build SOS problem for bisection.

opts = obj.opts;


% problem size
n = length(nlsos.x);
m = length(nlsos.g);

% get cone dimensions
if isfield(opts.Kx,'l'), Nl = opts.Kx.l; else, Nl = 0; end
if isfield(opts.Kx,'s'), Ns = opts.Kx.s; else, Ns = 0; end
if isfield(opts.Kc,'l'), Ml = opts.Kc.l; else, Ml = 0; end
if isfield(opts.Kc,'s'), Ms = opts.Kc.s; else, Ms = 0; end

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
p0 = sos.p;
        

% for each decision variable generate a parameter (which is the current solution)
% naming is not important since this is only for internal 
base_x = basis(sos.x);

x0    = casos.PS.sym('X0',base_x);
x1    = casos.PS.sym('X1',base_x);
lam_x = casos.PS.sym('lam_x',base_x);

% extend parameter vector
sos.p = [p0; x0];

% Taylor Approximation constraints and evaluate at parameterized
% solution
sos.g = linearize(sos.g,sos.x,x0);

% Taylor Approximation cost and evaluate at parameterized
% solution
sos.f = linearize(sos.f,sos.x,x0);

% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = struct('l',Nl,'s',Ns);
sosopt.Kc = struct('l',Ml,'s',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% see SOSOPTCOMMON#GET_MONOMIALS_IN for details
obj.monom_xl = monomials_in(obj.sossolver,2);
obj.monom_xs = monomials_in(obj.sossolver,4);
obj.monom_p  = basis(nlsos.p);
obj.monom_f  = basis(nlsos.f);
obj.monom_gl = monomials_in(obj.sossolver,5);
obj.monom_gs = monomials_in(obj.sossolver,7);


% setup line-search
d = casos.PS.sym('d');
  
ObjFun   = casos.Function('f',{sos.x}, {sos.f});
ConFun   = casos.Function('f',{sos.x}, {sos.g});

eta = 1e-15;


lam_g = casos.PS.sym('lam_g',obj.monom_gs);

xi  = x0*(1-d) + d*x1;
Psi = ObjFun(xi) - dot(lam_x,x1) - dot(lam_g, ConFun(xi));

% define SOS problem:   min_d Psi(d) s.t. 0 \leq d \leq 1 
sos_lineSearch = struct('x',d, ...
                        'f',Psi - eta*d, ...
                        'g',[], ... % 0 <= d <= 1  is set in call
                        'p',[x0; x1; lam_x; lam_g]);

       
% solve by relaxation to SDP
obj.lineSearch = casos.sossol('S','mosek',sos_lineSearch,...
                              struct('Kc',struct('s',0),'Kx',struct('l',1)));
       
         

end
