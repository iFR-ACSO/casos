function buildproblem(obj,sos)
% Build Alfonso problem from SOS dual formulation.
%
% Solves the primal affine SOS problem
%
%   min c(x) s.t. x = (xl,xc) and
%       Al(x) - bl in [lbg ubg] 
%       Ac(x) - bc in Kc
%       xl in [lbx ubx]
%       xc in Kx
%
% via its dual problem
%
%   max yb1(lbx) - yb2(ubx) + yl1(bl+lbg) - yl2(bl+ubg) + yc(bc) 
%       s.t. y = (yx,yl,yc) and
%       c - [yx1-yx2] - Al*(yl1-yl2) - Ac*(yc) in {0} x Kx*
%       yx, yl in R+ 
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
pts_lin_x = interpolation(Zvar_l);
pts_sos_x = interpolation(Zvar_s);
pts_lin_g = interpolation(Zcon_l);
pts_sos_g = interpolation(Zcon_s);

% compute Vandermonde matrix for decision variables
Vx_l = vandermonde2(Zvar_l,pts_lin_x);
Vx_s = vandermonde2(Zvar_s,pts_sos_x);
% and constraints
Vg_l = vandermonde2(Zcon_l,pts_lin_g);
Vg_s = vandermonde2(Zcon_s,pts_sos_g);
% combined Vandermonde matrices
Vx = blkdiag(Vx_l,Vx_s);
Vg = blkdiag(Vg_l,Vg_s);

% evaluate constraints on interpolation
Pcon = poly2interpol(sos.g,pts_sos_g,Zcon_s);

% number of linear variables / constraints
nnz_lin_x = nnz(Zvar_l);
nnz_lin_g = nnz(Zcon_l);
% number of sum-of-squares variables / constraints
nnz_sos_x = nnz(Zvar_s);
nnz_sos_g = nnz(Zcon_s);
% total number of variables / constraints
nnz_x = nnz_lin_x + nnz_sos_x;
nnz_g = nnz_lin_g + nnz_sos_g;

assert(length(Qvar) == nnz_x, 'Sum-of-squares decision variables must be in Gram form.')

% linearize cost and constraints w.r.t. monomial basis
c_monom = jacobian(Qobj,Qvar);
A_monom = jacobian(Pcon,Qvar); 
% output is already in interpolation basis
b_inter = -Pcon;
% convert input to interpolation basis
c_inter = c_monom/Vx;
A_inter = A_monom/Vx;

% linear bounds (monomial basis)
lbx = casadi.SX.sym('lbx',nnz_lin_x,1);
ubx = casadi.SX.sym('ubx',nnz_lin_x,1);
lbg = casadi.SX.sym('lbg',nnz_lin_g,1);
ubg = casadi.SX.sym('ubg',nnz_lin_g,1);
% convert bounds to interpolation basis
lbx_inter = Vx_l*lbx;
ubx_inter = Vx_l*ubx;
lbg_inter = Vg_l*lbg;
ubg_inter = Vg_l*ubg;

% separate linear and conic constraints
A = mat2cell(A_inter,[nnz_lin_g nnz_sos_g],nnz_x);
B = mat2cell(b_inter,[nnz_lin_g nnz_sos_g],1);

% build LMI for dual problem
%
%   min b(y) s.t. y = (yb,yl,yc) and
%       At(y) - c = [0; yx]
%       yx in Kx*
%       yc in Kc*
%
lmi.b = [                       % will be evaluated with x = 0
         -lbx_inter
          ubx_inter
         -B{1}-lbg_inter
          B{1}+ubg_inter
         -B{2}
];
lmi.At = [[
         -casadi.DM.eye(nnz_lin_x) +casadi.DM.eye(nnz_lin_x)
          casadi.DM(nnz_sos_x,2*nnz_lin_x)
] -A{1}' +A{1}' -A{2}'];
lmi.c  = -c_inter';
% half-basis for SOS dual cones
% [Px,~] = qr(vandermonde(Zvar_s0,pts_sos_x));
% [Pg,~] = qr(vandermonde(Zcon_s0,pts_sos_g));

% remove decision variables
fopt = struct('allow_free',true);
lmi = call(casadi.Function('lmi',{Qvar Qpar},struct2cell(lmi),{'x' 'p'},fieldnames(lmi),fopt),struct('x',0,'p',Qpar));

% Alfonso solves the primal problem
%
%   min c(x) s.t. x = (yb,yl,yc,yx) and
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
% set cones
cone = [
    % linear cone (corresponds to linear polynomial constraints)
    {struct('type','l','dim',2*nnz_lin_x+2*nnz_lin_g);}
    % LMI cones 1 (corresponds to SOS constraints)
    arrayfun(@(i) get_dualcone(Zcon_s,Zcon_s0,pts_sos_g,i), 1:Ms, 'UniformOutput', false)'
    % LMI cone 2 (corresponds to SOS variables)
    arrayfun(@(i) get_dualcone(Zvar_s,Zvar_s0,pts_sos_x,i), 1:Ns, 'UniformOutput', false)'
];
% cone{2} = struct('type','grk1lmi','dim',nnz_sos_g,'nonneg',false,'ext',false);
% cone{2}.Vs = {full(Pg')};
% cone{2}.ws = ones(nnz_sos_g,1);
% cone{2}.x0 = ones(nnz_sos_g,1);
% cone{3} = struct('type','grk1lmi','dim',nnz_sos_x,'nonneg',false,'ext',false);
% cone{3}.Vs = {full(Px')};
% cone{3}.ws = ones(nnz_sos_x,1);
% cone{3}.x0 = ones(nnz_sos_x,1);

% save problem data and cones
obj.fhan = casadi.Function('f', ...
                            {Qpar lbx ubx lbg ubg}, ...
                            struct2cell(prob), ...
                            {'p' 'lbx' 'ubx' 'lbg' 'ubg'}, ...
                            fieldnames(prob) ...
);
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
lmisol.x = casadi.SX.sym('lmi_x',nnz_x+nnz_g+nnz_lin_x+nnz_lin_g);
lmisol.s = casadi.SX.sym('lmi_s',nnz_x+nnz_g+nnz_lin_x+nnz_lin_g);
lmisol.y = casadi.SX.sym('lmi_y',nnz_x);
% lower constraint bounds
lmisol.lbg = lbg;

% separate primal and dual slack variables (of LMI)
Lam = mat2cell(lmisol.x,[nnz_lin_x nnz_lin_x nnz_lin_g nnz_lin_g nnz_sos_g nnz_sos_x]);
Con = mat2cell(lmisol.s,[nnz_lin_x nnz_lin_x nnz_lin_g nnz_lin_g nnz_sos_g nnz_sos_x]);

% coordinates of SOS solution
sossol.x = Vx\lmisol.y;
sossol.f = -lmisol.f;
sossol.g = Vg\[Con{3}+lbg_inter; Con{5}];
sossol.lam_x = -Vx'*[Lam{1}-Lam{2}; Lam{6}];
sossol.lam_g = -Vg'*[Lam{3}-Lam{4}; Lam{5}];

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

function cone = get_dualcone(Z_gram,Z_half,pts,i)
% Return dual cone structure.

    % number of terms (nonzeros) in gram basis
    nT = Z_gram.get_nterm(i);
    % get Vandermonde matrix
    V = vandermonde2(Z_half,pts(:,1:nT),i);
    % orthonormalize
    [P,~] = qr(V);

    % set LMI cone
    cone = struct('type','grk1lmi','dim',nT,'nonneg',false,'ext',false);
    cone.Vs = {full(P')};
    cone.ws = ones(nT,1);
    cone.x0 = ones(nT,1);
end
