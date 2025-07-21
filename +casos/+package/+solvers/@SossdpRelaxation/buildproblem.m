function buildproblem(obj,solver,sos)
% Build SDP problem from SOS relaxation.

opts   = obj.opts;

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
[Zvar_s,Ksdp_x_s,Zvar_s0,Mp_x,Md_x] = grambasis(sparsity(sos.x),Is,opts.newton_solver);
% obtain Gram basis for sum-of-squares constraints
[Zcon_s,Ksdp_g_s,Zcon_s0,Mp_g,Md_g] = grambasis(sparsity(sos.g),Js,opts.newton_solver);

% matrix decision variables for variables
Qvar_G = casadi.SX.sym('P',sum(Ksdp_x_s.^2),1);
% matrix decision variables for constraints
Qcon_G = casadi.SX.sym('Q',sum(Ksdp_g_s.^2),1);

% linear decision variables
[Qvar_l,Zvar_l] = poly2basis(sos.x,~Is); % TODO: internal
% linear constraints
Zcon_l = basis(sos.g,~Js);
% handle parameters
[Qpar,Zpar] = poly2basis(sos.p);
% get cost function
[Qobj,Zobj] = poly2basis(sos.f);

% get combined variables / constraints
[Qvar,Zvar] = poly2basis(sos.x);
% [Qcon,Zcon] = poly2basis(sos.g,[Zcon_l; Zcon_s]);
Zcon = [Zcon_l; Zcon_s];

% number of linear variables / constraints
nnz_lin_x = nnz(Zvar_l);
nnz_lin_g = nnz(Zcon_l);
% number of sum-of-squares variables / constraints
nnz_sos_x = nnz(Zvar_s);
nnz_sos_g = nnz(Zcon_s);
% number of Gram variables / constraints
nnz_gram_x = sum(Ksdp_x_s.^2);
nnz_gram_g = sum(Ksdp_g_s.^2);

assert(length(Qvar) == (nnz_lin_x + nnz_sos_x), 'Sum-of-squares decision varibles must be in Gram form.')

if strcmp(opts.basis_type, 'interpolation')
    %% Relaxation to interpolation basis
    % generate interpolation basis
    % pts_lin_x = interpolation(Zvar_l);
    pts_sos_x = interpolation(Zvar_s);
    % pts_lin_g = interpolation(Zcon_l);
    pts_sos_g = interpolation(Zcon_s);
    
    % compute Vandermonde matrix for decision variables
    Vx_l = speye(nnz_lin_x);
    % Vx_l = vandermonde2(Zvar_l,pts_lin_x);
    Vx_s = vandermonde2(Zvar_s,pts_sos_x);
    % and constraints
    Vg_l = speye(nnz_lin_g);
    % Vg_l = vandermonde2(Zcon_l,pts_lin_g);
    Vg_s = vandermonde2(Zcon_s,pts_sos_g);
    % combined Vandermonde matrices
    % Vx = blkdiag(Vx_l,Vx_s);
    % Vg = blkdiag(Vg_l,Vg_s);

    % evaluate constraints on interpolation
    assert(nnz_lin_g == 0, 'Not implemented.') % TODO
    Qcon = poly2interpol(sos.g,pts_sos_g,Zcon_s);

    % compute map from SDP to interpolation for each element
    WW_x = arrayfun(@(i) interp_map(Zvar_s,Zvar_s0,pts_sos_x,i), 1:Ns, 'UniformOutput', false);
    WW_g = arrayfun(@(i) interp_map(Zcon_s,Zcon_s0,pts_sos_g,i), 1:Ms, 'UniformOutput', false);
    % combine
    Mp_x = blkdiag(WW_x{:});
    Mp_g = blkdiag(WW_g{:});

    % compute dual maps as left-inverse
    Md_x = (Mp_x*Mp_x')\Mp_x;
    Md_g = (Mp_g*Mp_g')\Mp_g;

else
    %% Relaxation to monomial basis
    % get constraints on monomial basis
    Qcon = poly2basis(sos.g,Zcon);

    % weighting matrices in lieu of Vandermonde
    Vx_s = speye(nnz_sos_x);
    Vg_s = speye(nnz_sos_g);
end

% matrix decision variables
Qvar_sdp = [Qvar_l; Qvar_G];

% replace sum-of-squares decision variables
% and pre-compute derivatives
% gradient and hessian of objective
sosprob_f = casadi.Function('sos_f',{Qvar},{hessian(Qobj,Qvar) jacobian(Qobj,Qvar) Qobj},struct('allow_free',true)); %hessian_old(sosprob,0,0);
% jacobian of constraints
sosprob_g = casadi.Function('sos_g',{Qvar},{jacobian(Qcon,Qvar) Qcon},struct('allow_free',true)); %jacobian_old(sosprob,0,1);

% map sum-of-squares decision variables to matrix variables
map_x = blkdiag(speye(nnz_lin_x), Vx_s\Mp_x);
Qvar_mapped = map_x*Qvar_sdp;
% [sdp_f,sdp_g] = sosprob(Qvar_mapped);
[sdp_Hf,sdp_Jf,sdp_f] = sosprob_f(Qvar_mapped);
[sdp_Jg,sdp_g] = sosprob_g(Qvar_mapped);

% map sum-of-squares constraints to matrix variables
map_g = blkdiag(sparse(nnz_lin_g,0), Mp_g);

% build SDP problem
sdp.x = [Qvar_sdp; Qcon_G];
sdp.f = sdp_f;
sdp.g = sdp_g - map_g*Qcon_G;
sdp.p = Qpar;
% store derivatives
sdp.derivatives.Hf = blockcat(map_x'*sdp_Hf*map_x, ...
                              sparse(nnz_lin_x+nnz_gram_x,nnz_gram_g), ...
                              sparse(nnz_gram_g,nnz_lin_x+nnz_gram_x), ...
                              sparse(nnz_gram_g,nnz_gram_g));
sdp.derivatives.Jf = horzcat(sdp_Jf*map_x, sparse(1,nnz_gram_g));
sdp.derivatives.Jg = horzcat(sdp_Jg*map_x, -map_g);
% SDP options
sdpopt = opts.sdpsol_options;
sdpopt.Kx = struct('lin', nnz_lin_x, 'psd', [Ksdp_x_s; Ksdp_g_s]);
sdpopt.Kc = struct('lin', nnz_lin_g + nnz_sos_g);

% initialize SDP solver
obj.sdpsolver = casos.package.solvers.SdpsolInternal('SDP',solver,sdp,sdpopt);

% store basis
obj.sparsity_x  = Zvar;
obj.sparsity_xl = Zvar_l;
obj.sparsity_xs = Zvar_s;
obj.sparsity_p  = Zpar;
obj.sparsity_f  = Zobj;
obj.sparsity_g  = Zcon;
obj.sparsity_gl = Zcon_l;
obj.sparsity_gs = Zcon_s;

% map SDP solution to SOS solution
% symbolic SDP solution
sdpsol.x = casadi.SX.sym('sol_x',size(sdp.x));
sdpsol.f = casadi.SX.sym('sol_f');
sdpsol.g = casadi.SX.sym('sol_g',size(sdp.g));
sdpsol.lam_x = casadi.SX.sym('sol_lam_x',size(sdp.x));
sdpsol.lam_g = casadi.SX.sym('sol_lam_g',size(sdp.g));
sdpsol.lam_p = casadi.SX.sym('sol_lam_p',size(sdp.p));

% coordinates of SOS solution
sossol.x = blkdiag(speye(nnz_lin_x), Vx_s\Mp_x, sparse(0,nnz_gram_g))*sdpsol.x;
sossol.f = sdpsol.f;
sossol.g = [
    blkdiag(speye(nnz_lin_g), sparse(0,nnz_sos_g))*sdpsol.g
    blkdiag(sparse(0,nnz_lin_x+nnz_gram_x), Vg_s\Mp_g)*sdpsol.x
];
sossol.lam_x = blkdiag(speye(nnz_lin_x), Vx_s'*Md_x, sparse(0,nnz_gram_g))*sdpsol.lam_x;
sossol.lam_g = [
    blkdiag(speye(nnz_lin_g), sparse(0,nnz_sos_g))*sdpsol.lam_g
    blkdiag(sparse(0,nnz_lin_x+nnz_gram_x),  Vg_s'*Md_g)*sdpsol.lam_x
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

function W = interp_map(Z_gram,Z_half,pts,i)
% Return map from SDP to interpolation basis.
%
% This computes the adjoint of the operator
%
%   y -> P'*diag(y)*P
%
% where P is the ortonormalization of the Vandermonde matrix for the Gram
% half-basis.

    % number of terms (nonzeros) in gram basis
    nT = Z_gram.get_nterm(i);
    % get Vandermonde matrix
    V = vandermonde2(Z_half,pts(:,1:nT),i);
    % orthonormalize
    [P,~] = qr(V);

    % Kronecker product
    VV = kron(P,P);
    % indices into terms of Gram basis
    k = 1:nT;
    % index into Kronecker product
    idx = (nT+1)*(k-1)+1;
    W = VV(idx,:);
end
