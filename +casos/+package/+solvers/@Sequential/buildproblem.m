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
        
% parameterize decision variables
x0    = [];
x1    = [];
lam_x = [];
lam_g = [];

% for each decision variable generate a parameter (which is the current solution)
% naming is not important since this is only for internal 
for k = 1:length(sos.x)

     base = basis(sos.x(k));
     x0 = [x0; 
           casos.PS.sym(['X0' num2str(k)],base)];

     x1 = [x1; 
           casos.PS.sym(['X1' num2str(k)],base)];

     lam_x = [lam_x; 
              casos.PS.sym(['lam_x' num2str(k)],base)];
end


for k = 1:length(sos.g)
     lam_g = [lam_g; 
              casos.PS.sym(['lam_g' num2str(k)],basis(sos.g(k)))];
end

        
d = casos.PS.sym('d');
  
ObjFun   = casos.Function('f',{sos.x}, {sos.f});
ConFun   = casos.Function('f',{sos.x}, {sos.g});

eta = 1e-15;


% setup line search solver

% define SOS problem:   min_d Psi(d) s.t. 0 \leq d \leq 1 
sos_lineSearch = struct('x',d, ...
                        'f',ObjFun((1-d)*x0 + d*x1) - dot(lam_x,x0) - dot(lam_g, ConFun((1-d)*x0 + d*x1)) - eta*d.^2, ...
                        'g',[], ...
                        'p',[x0; x1; lam_x; lam_g]);

       
% solve by relaxation to SDP
obj.lineSearch = casos.sossol('S','mosek',sos_lineSearch,...
                              struct('Kc',struct('s',0),'Kx',struct('l',1)));
       
         
% Taylor Approximation constraints and evaluate at parameterized
% solution
sos.g = linearize(sos.g,sos.x,x0);

% Taylor Approximation cost and evaluate at parameterized
% solution
sos.f = linearize(sos.f,sos.x,x0);
        
% extend parameter vector
sos.p = [p0; x0];

% SOS options
sosopt = opts.sossol_options;
sosopt.Kx = struct('l',Nl,'s',Ns);
sosopt.Kc = struct('l',Ml,'s',Ms);
sosopt.error_on_fail = false;

% initialize convex SOS solver
obj.sossolver = casos.package.solvers.SossolInternal('SOS',opts.sossol,sos,sosopt);

% store basis
obj.monom_xl = basis(nlsos.x,~Is);
obj.monom_xs = basis(nlsos.x, Is);
obj.monom_p  = basis(nlsos.p);
obj.monom_f  = basis(nlsos.f);
obj.monom_gl = basis(nlsos.g,~Js);
obj.monom_gs = basis(nlsos.g, Js);

end
