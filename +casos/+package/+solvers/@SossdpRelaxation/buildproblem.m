function buildproblem(obj,solver,sos)
% Build SDP problem from SOS relaxation.

opts = obj.opts;

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

% obtain Gram basis for decision variables
[Zgram_x,Ksdp_x_s,~,Mp_x,Md_x] = grambasis(sparsity(sos.x),Is);
% obtain Gram basis for sum-of-squares constraints
[Zgram_g,Ksdp_g_s,~,Mp_g,Md_g] = grambasis(sparsity(sos.g),Js);

% matrix decision variables for variables
Qgram_x = casadi.SX.sym('P',sum(Ksdp_x_s.^2),1);
% matrix decision variables for constraints
Qgram_g = casadi.SX.sym('Q',sum(Ksdp_g_s.^2),1);

% handle decision variables
[Qlin_x,Zlin_x] = poly2basis(sos.x(~Is)); % TODO: internal
[Qsos_x,Zsos_x] = poly2basis(sos.x(Is),Zgram_x);
% handle constraints
[Qlin_g,Zlin_g] = poly2basis(sos.g(~Js)); % TODO: internal
[Qsos_g,Zsos_g] = poly2basis(sos.g(Js),Zgram_g);
% handle parameter
[Qlin_p,Zlin_p] = poly2basis(sos.p);

% get cost function
[Qlin_f,Zlin_f] = poly2basis(sos.f);

% number of linear variables / constraints
nnz_lin_x = nnz(Zlin_x);
nnz_lin_g = nnz(Zlin_g);
% number of sum of squares variables / constraints
nnz_sos_x = nnz(Zsos_x);
nnz_sos_g = nnz(Zsos_g);
% number of Gram variables / constraints
nnz_gram_x = sum(Ksdp_x_s.^2);
nnz_gram_g = sum(Ksdp_g_s.^2);

% replace sum-of-squares decision variables
Df_x = jacobian(Qlin_f, [Qlin_x;Qsos_x]);
Dg_x = jacobian([Qlin_g;Qsos_g],[Qlin_x;Qsos_x]);
% map sum-of-squares decision variables to matrix variables
Df_x_mapped = Df_x*blkdiag(speye(nnz_lin_x), Mp_x);
Dg_x_mapped = Dg_x*blkdiag(speye(nnz_lin_x), Mp_x);

% constant objective / constraint
f0 = mtaylor(Qlin_f, [Qlin_x;Qsos_x], 0, 0);
g0 = mtaylor([Qlin_g;Qsos_g], [Qlin_x;Qsos_x], 0, 0);

% build SDP problem
sdp.x = [Qlin_x; Qgram_x; Qgram_g];
sdp.f = f0 + Df_x_mapped*[Qlin_x; Qgram_x];
sdp.g = g0 + Dg_x_mapped*[Qlin_x; Qgram_x] - Mp_g*Qgram_g;
% sdp.f = Qlin_f;
% sdp.g = Qdiff_g;
sdp.p = Qlin_p;
% SDP options
sdpopt = opts.sdpsol_options;
sdpopt.Kx = struct('lin', nnz_lin_x, 'psd', [Ksdp_x_s; Ksdp_g_s]);
sdpopt.Kc = struct('lin', nnz_lin_g + nnz_sos_g);

% initialize SDP solver
obj.sdpsolver = casos.package.solvers.SdpsolInternal('SDP',solver,sdp,sdpopt);

% store basis
obj.sparsity_xl = Zlin_x;
obj.sparsity_xs = Zsos_x;
obj.sparsity_p  = Zlin_p;
obj.sparsity_f  = Zlin_f;
obj.sparsity_gl = Zlin_g;
obj.sparsity_gs = Zsos_g;

% map SDP solution to SOS solution
% symbolic SDP solution
sdpsol.x = casadi.SX.sym('sol_x',size(sdp.x));
sdpsol.f = casadi.SX.sym('sol_f');
sdpsol.g = casadi.SX.sym('sol_g',size(sdp.g));
sdpsol.lam_x = casadi.SX.sym('sol_lam_x',size(sdp.x));
sdpsol.lam_g = casadi.SX.sym('sol_lam_g',size(sdp.g));
sdpsol.lam_p = casadi.SX.sym('sol_lam_p',size(sdp.p));

% coordinates of SOS solution
sossol.x = blkdiag(speye(nnz_lin_x), Mp_x, sparse(0,nnz_gram_g))*sdpsol.x;
sossol.f = sdpsol.f;
sossol.g = [
    blkdiag(speye(nnz_lin_g), sparse(0,nnz_sos_g))*sdpsol.g
    blkdiag(sparse(0,nnz_lin_x+nnz_gram_x),  Mp_g)*sdpsol.x
];
sossol.lam_x = blkdiag(speye(nnz_lin_x), Md_x, sparse(0,nnz_gram_g))*sdpsol.lam_x;
sossol.lam_g = [
    blkdiag(speye(nnz_lin_g), sparse(0,nnz_sos_g))*sdpsol.lam_g
    blkdiag(sparse(0,nnz_lin_x+nnz_gram_x),  Md_g)*sdpsol.lam_x
];

% options for Casadi functions
% NOTE: As per Casadi v3.6.5, functions' input and output names 
% must be mutually exclusive unless explicity permitted.
fopt = struct('allow_duplicate_io_names',true);

% output function
obj.gram2sos = casadi.Function('L', ...
                struct2cell(sdpsol), ...
                struct2cell(sossol), ...
                fieldnames(sdpsol), ...
                fieldnames(sossol), ...
                fopt ...
);

end
