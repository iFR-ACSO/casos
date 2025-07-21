function buildproblem(obj,sos)
% Build LMI relaxation from SOS dual formulation.
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
%       s.t. y = (yb,yl,yc) and
%       c - [yb1-yb2; 0] - Al*(yl1-yl2) - Ac*(yc) in {0} x Kx*
%       yb, yl in R+ 
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

% separate linear and conic constraints + variables
A = mat2cell(A_inter,[nnz_lin_g nnz_sos_g],[nnz_lin_x nnz_sos_x]);
B = mat2cell(b_inter,[nnz_lin_g nnz_sos_g],1);
C = mat2cell(c_inter,1,[nnz_lin_x nnz_sos_x]);

% build LMI for dual problem
%
%   min b(y) s.t. y = (yb,yl,yc) and
%       c - At(y) = [0; yx]
%       yx in Kx*
%       yc in Kc*
%
yb1 = casadi.SX.sym('yb1',nnz_lin_x,1);
yb2 = casadi.SX.sym('yb2',nnz_lin_x,1);
yl1 = casadi.SX.sym('yl1',nnz_lin_g,1);
yl2 = casadi.SX.sym('yl2',nnz_lin_g,1);
yx  = casadi.SX.sym('yx', nnz_sos_x,1);
yc  = casadi.SX.sym('yc', nnz_sos_g,1);

% constraints (corresponds to linear variables)
gl = C{1}' - (yb1-yb2) - A{1,1}'*(yl1-yl2) - A{2,1}'*yc;
% constraints (corresponds to SOS variables)
gs = C{2}' - yx        - A{1,2}'*(yl1-yl2) - A{2,2}'*yc;

% build LMI
%
%   min f(x) s.t. x = (yb,yl,yc,yx) and
%       g(x) = 0
%       Pc'(diag yc)Pc is p.s.d.
%       Px'(diag yx)Px is p.s.d.
lmi.f = -lbx_inter'*yb1 + ubx_inter'*yb2 - (B{1}+lbg_inter)'*yl1 + (B{1}+ubg_inter)'*yl2 - B{2}'*yc;

% number of terms in gram basis
nterm_gram_x = get_nterm(Zvar_s);
nterm_gram_g = get_nterm(Zcon_s);

% separate decision variables
Yx = mat2cell(yx,nterm_gram_x,1);
Yc = mat2cell(yc,nterm_gram_g,1);
% LMI representation of dual cones
G = [
    % LMI cones 1 (corresponds to SOS constraints)
    arrayfun(@(i) get_dualcone(Yc{i},Zcon_s0,nterm_gram_g(i),pts_sos_g,i), 1:Ms, 'UniformOutput', false)'
    % LMI cone 2 (corresponds to SOS variables)
    arrayfun(@(i) get_dualcone(Yx{i},Zvar_s0,nterm_gram_x(i),pts_sos_x,i), 1:Ns, 'UniformOutput', false)'
];
% combine constraints
lmi.g = vertcat(gl, gs, G{:});
% decision variables
lmi.x = [yb1;yb2;yl1;yl2;yc;yx];
% parameters
lmi.p = [Qpar;Qvar;lbx;ubx;lbg;ubg];

% TODO: precompute derivatives
VV = jacobian(vertcat(G{:}),[yc;yx]);

% dimensions of dual cones
dim_cone_x = get_nterm(Zvar_s0);
dim_cone_g = get_nterm(Zcon_s0);

% options to LMI solver
sdpopt = opts.sdpsol_options;
sdpopt.Kx = struct('lin', 2*nnz_lin_x+2*nnz_lin_g+nnz_sos_g+nnz_sos_x);
sdpopt.Kc = struct('lin', nnz_lin_x+nnz_sos_x, 'psd', full([dim_cone_g; dim_cone_x]));

% initialize LMI solver
obj.sdpsolver = casos.package.solvers.SdpsolInternal('LMI',opts.sdpsol,lmi,sdpopt);

% options for Casadi functions
% NOTE: As per Casadi v3.6.5, functions' input and output names 
% must be mutually exclusive unless explicity permitted.
fopt = struct('allow_duplicate_io_names',true);

% arguments to LMI solver
obj.lmiargs = casadi.Function('f', {Qpar lbx ubx lbg ubg}, { ...
                [Qpar;zeros(size(Qvar));lbx;ubx;lbg;ubg] % parameters to LMI
                [zeros(2*nnz_lin_x+2*nnz_lin_g,1); -inf(nnz_sos_g+nnz_sos_x,1)] % lower bound LMI variables
                inf(2*nnz_lin_x+2*nnz_lin_g+nnz_sos_g+nnz_sos_x,1)              % upper bound LMI variables
                0       % lower bound LMI constraints
                0       % upper bound LMI constraints
},{'p' 'lbx' 'ubx' 'lbg' 'ubg'},{'p' 'lbx' 'ubx' 'lbg' 'ubg'},fopt);

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
lmisol.lam_x = casadi.SX.sym('lmi_lam_x',nnz_x+nnz_g+nnz_lin_x+nnz_lin_g);
lmisol.lam_g = casadi.SX.sym('lmi_lam_g',nnz_lin_x+nnz_sos_x+sum([dim_cone_g; dim_cone_x].^2));
% lmisol.s = casadi.SX.sym('lmi_s',nnz_x+nnz_g+nnz_lin_x+nnz_lin_g);
% lmisol.y = casadi.SX.sym('lmi_y',nnz_x);
lmisol.g = casadi.SX.sym('lmi_g',size(lmi.g));
lmisol.lam_p = casadi.SX.sym('lmi_lam_p',size(lmi.p));
% lower constraint bounds
lmisol.lbg = lbg;

lmisol_y = -lmisol.lam_g(1:nnz_x);
lmisol_s = [-lmisol.lam_x(1:(2*nnz_lin_x+2*nnz_lin_g)); -VV'*lmisol.lam_g(nnz_x+1:end)];

% separate primal and dual slack variables (of LMI)
Lam = mat2cell(lmisol.x,[nnz_lin_x nnz_lin_x nnz_lin_g nnz_lin_g nnz_sos_g nnz_sos_x]);
Con = mat2cell(lmisol_s,[nnz_lin_x nnz_lin_x nnz_lin_g nnz_lin_g nnz_sos_g nnz_sos_x]);

% coordinates of SOS solution
sossol.x = Vx\lmisol_y;
sossol.f = -lmisol.f;
sossol.g = Vg\[Con{3}+lbg_inter; Con{5}];
sossol.lam_x = -Vx'*[Lam{1}-Lam{2}; Lam{6}];
sossol.lam_g = -Vg'*[Lam{3}-Lam{4}; Lam{5}];

% options for Casadi functions
% NOTE: As per Casadi v3.6.5, functions' input and output names 
% must be mutually exclusive unless explicity permitted.
fopt = struct('allow_duplicate_io_names',true);

% output function
obj.lmi2sos = casadi.Function('g', ...
                struct2cell(lmisol), ...
                struct2cell(sossol), ...
                fieldnames(lmisol), ...
                fieldnames(sossol), ...
                fopt ...
);

end

function cone = get_dualcone(y,Z_half,nT,pts,i)
% Return dual cone structure.

    % get Vandermonde matrix
    V = vandermonde2(Z_half,pts(:,1:nT),i);
    % orthonormalize
    [P,~] = qr(V);

    % set LMI cone
    cone = reshape(P'*diag(y)*P,size(P,2).^2,1);
end
