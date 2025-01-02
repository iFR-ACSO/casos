function buildproblem(obj,sos)
% Build Alfonso problem from SOS dual formulation.
%
% Solves the primal affine SOS problem
%
%   min c(x) s.t. x = (xl,xc) and
%       A(x) + b in {0} x Kc
%       xl in R[.]
%       xc in Kx
%
% via its dual problem
%
%   max y(-b) s.t. y = (yl,yc) and
%       c - A*(y) in {0} x Kx*
%       yl in R[.]*
%       yc in Kc*
%
% where the primal decision variables and the output of the linear operator
% are given in the interpolation basis.

opts   = obj.opts;
newton = opts.newton_simplify;

% problem size
n = length(sos.x);
m = length(sos.g);

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

% linear decision variables
Zvar_l = basis(sos.x,~Is);
% linear constraints
Zcon_l = basis(sos.g,~Js);
% handle parameters
[Qpar,Zpar] = poly2basis(sos.p);
% get cost function
[Qobj,Zobj] = poly2basis(sos.f);

% obtain Gram basis for decision variables
[Zvar_s,~,Zvar_s0] = grambasis(sparsity(sos.x),Is,newton);
% obtain Gram basis for sum-of-squares constraints
[Zcon_s,~,Zcon_s0] = grambasis(sparsity(sos.g),Js,newton);

% get combined variables / constraints
[Qvar,Zvar] = poly2basis(sos.x);
Zcon = [Zcon_l; Zcon_s];

% generate interpolation basis
pts_x = interpolation(Zvar);
pts_g = interpolation(Zcon);

% compute Vandermonde matrix for decision variables
Vx = vandermonde(Zvar,pts_x);
% and constraints
Vg = vandermonde(Zcon,pts_g);

% evaluate constraints on interpolation
Pcon = poly2interpol(sos.g,pts_g);

% number of linear variables / constraints
nnz_lin_x = nnz(Zvar_l);
nnz_lin_g = nnz(Zcon_l);
% number of sum-of-squares variables / constraints
nnz_sos_x = nnz(Zvar_s);
nnz_sos_g = nnz(Zcon_s);
% number of Gram variables / constraints
% nnz_gram_x = sum(Ksdp_x_s.^2);
% nnz_gram_g = sum(Ksdp_g_s.^2);

assert(length(Qvar) == (nnz_lin_x + nnz_sos_x), 'Sum-of-squares decision variables must be in Gram form.')

% linearize cost and constraints w.r.t. monomial basis
c_monom = jacobian(Qobj,Qvar);
A_monom = jacobian(Pcon,Qvar); % output is already in interpolation basis
% convert input to interpolation basis
c_inter = c_monom/Vx;
A_inter = A_monom/Vx;

% build LMI for dual problem
%
%   min b(y) s.t. y = (yl,yc) and
%       At(y) - c = [0; yx]
%       yx in Kx*
%       yc in Kc*
%
lmi.b = Pcon;   % will be evaluated with x = 0
lmi.At = -A_inter';
lmi.c  = -c_inter';
% half-basis for SOS dual cones
[Px,~] = qr(vandermonde(Zvar_s0,pts_x));
[Pg,~] = qr(vandermonde(Zcon_s0,pts_g));

% remove decision variables
lmi = call(casadi.Function('lmi',{Qvar Qpar},struct2cell(lmi),{'x' 'p'},fieldnames(lmi)),struct('x',0,'p',Qpar));

% Alfonso solves the primal problem
%
%   min c(x) s.t. x = (yl,yc,yx) and
%       A(x) = b
%       Pc'(diag yc)Pc is p.s.d.
%       Px'(diag yx)Px is p.s.d.
%
% with dual problem
%
%   max b(y) s.t. 
%       c - At(y) = s
%       s in dual cone
% 
prob.A = [lmi.At [casadi.DM(nnz_lin_x,nnz_sos_x); -casadi.DM.eye(nnz_sos_x)]];
prob.b = lmi.c;
prob.c = [lmi.b; casadi.DM(nnz_sos_x,1)];
% linear cone (corresponds to linear polynomial constraints)
cone{1} = struct('type','free','dim',nnz_lin_g);
% LMI cone 1 (corresponds to SOS constraints)
cone{2} = struct('type','grk1lmi','dim',nnz_sos_g,'nonneg',false,'ext',false);
cone{2}.Vs = {full(Pg')};
cone{2}.ws = ones(nnz_sos_g,1);
cone{2}.x0 = ones(nnz_sos_g,1);
% LMI cone 2 (corresponds to SOS variables)
cone{3} = struct('type','grk1lmi','dim',nnz_sos_x,'nonneg',false,'ext',false);
cone{3}.Vs = {full(Px')};
cone{3}.ws = ones(nnz_sos_x,1);
cone{3}.x0 = ones(nnz_sos_x,1);

% save problem data and cones
obj.fhan = casadi.Function('f',{Qpar},struct2cell(prob),{'p'},fieldnames(prob));
obj.cone = cone;

% store basis
obj.sparsity_x  = Zvar;
obj.sparsity_xl = Zvar_l;
obj.sparsity_xs = Zvar_s;
obj.sparsity_p  = Zpar;
obj.sparsity_f  = Zobj;
obj.sparsity_g  = Zcon;
obj.sparsity_gl = Zcon_l;
obj.sparsity_gs = Zcon_s;

% map LMI solution to SOS solution
% symbolic LMI solution
lmisol.f = casadi.SX.sym('lmi_f');
lmisol.x = casadi.SX.sym('lmi_x',nnz_lin_g+nnz_sos_g+nnz_sos_x);
lmisol.s = casadi.SX.sym('lmi_s',nnz_lin_g+nnz_sos_g+nnz_sos_x);
lmisol.y = casadi.SX.sym('lmi_y',nnz_lin_x+nnz_sos_x);

% separate primal and dual slack variables (of LMI)
Lam = mat2cell(lmisol.x,[nnz_lin_g+nnz_sos_g nnz_sos_x]);
Con = mat2cell(lmisol.s,[nnz_lin_g+nnz_sos_g nnz_sos_x]);

% coordinates of SOS solution
sossol.x = Vx\lmisol.y;
sossol.f = -lmisol.f;
sossol.g = Vg\Con{1};
sossol.lam_x = -Vx'*[zeros(nnz_lin_x,1); Lam{2}];
sossol.lam_g = -Vg'*Lam{1};

% options for Casadi functions
% NOTE: As per Casadi v3.6.5, functions' input and output names 
% must be mutually exclusive unless explicity permitted.
fopt = struct('allow_duplicate_io_names',true);

% output function
obj.ghan = casadi.Function('g', ...
            struct2cell(lmisol), ...
            struct2cell(sossol), ...
            fieldnames(lmisol), ...
            fieldnames(sossol), ...
            fopt ...
);

end
