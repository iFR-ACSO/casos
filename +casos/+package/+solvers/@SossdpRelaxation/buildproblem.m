function buildproblem(obj,solver,sos)
% Build SDP problem from SOS relaxation.

% extract options
opts   = obj.opts;

% problem size
n = length(sos.x);
m = length(sos.g);

% get cone dimensions
% get cone dimensions for the decision variables 
Nl   = get_dimension(obj.get_cones,opts.Kx,'lin');
Ns   = get_dimension(obj.get_cones,opts.Kx,'sos');
Nds  = get_dimension(obj.get_cones,opts.Kx,'dsos');
Nsds = get_dimension(obj.get_cones,opts.Kx,'sdsos');

% get cone dimensions for the constraints
Ml   = get_dimension(obj.get_cones,opts.Kc,'lin');
Ms   = get_dimension(obj.get_cones,opts.Kc,'sos');
Mds  = get_dimension(obj.get_cones,opts.Kc,'dsos');
Msds = get_dimension(obj.get_cones,opts.Kc,'sdsos');

% validate dimensions
assert(n == (Nl + Ns + Nds + Nsds), 'Dimension of Kx must be equal to number of variables (%d).', n);
assert(m == (Ml + Ms + Mds + Msds), 'Dimension of Kc must be equal to number of constraints (%d).', m)

% select sum-of-squares variables and constraints
Is = [false(Nl,1); true(Ns,1); true(Nds,1); true(Nsds,1)];
Js = [false(Ml,1); true(Ms,1); true(Mds,1); true(Msds,1)];

% obtain Gram basis for decision variables
[Zvar_s,Ksdp_x_s,~,Mp_x,Md_x] = grambasis(sparsity(sos.x),Is,opts.newton_solver);

% obtain Gram basis for sum-of-squares constraints
[Zcon_s,Ksdp_g_s,~,Mp_g,Md_g] = grambasis(sparsity(sos.g),Js,opts.newton_solver);

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
[Qcon,Zcon] = poly2basis(sos.g,[Zcon_l; Zcon_s]);

% number of linear variables / constraints
nnz_lin_x = nnz(Zvar_l);
nnz_lin_g = nnz(Zcon_l);

% number of sum-of-squares variables / constraints
nnz_sos_x = nnz(Zvar_s);
nnz_sos_g = nnz(Zcon_s);

% number of Gram variables / constraints
nnz_gram_x = sum(Ksdp_x_s.^2);
nnz_gram_g = sum(Ksdp_g_s.^2);

% validate sum-of-squares variable count
assert(length(Qvar) == (nnz_lin_x + nnz_sos_x), 'Sum-of-squares decision varibles must be in Gram form.')

% matrix decision variables
Qvar_sdp = [Qvar_l; Qvar_G];

% Create reordering map if needed (for DSOS/SDSOS)
if Nds + Nsds + Mds + Msds >= 0
    [permMat, ~, ~] = obj.process_reorder(Qvar_G, Qcon_G,     ...
                                          Ksdp_x_s, Ksdp_g_s, ...
                                          Qvar_l,             ...
                                          Ns, Ms, Nds, Mds);
else
    permMat = speye(length(Qvar_sdp) + length(Qcon_G));
end

% replace sum-of-squares decision variables
% and pre-compute derivatives
% gradient and hessian of objective
if ~isfield(sos,'derivatives')
    % compute derivatives if not  pre-computed
    Hf = hessian(Qobj,Qvar);
    Df = jacobian(Qobj,Qvar);
    Dg = jacobian(Qcon,Qvar);
else
    % use pre-computed derivatives (undocumented)
    Hf = op2basis(sos.derivatives.Hf,Zvar,Zvar);
    Df = op2basis(sos.derivatives.Df,Zvar,Zobj);
    Dg = op2basis(sos.derivatives.Dg,Zvar,Zcon);
end

% create CasADi functions for objective and constraints
sosprob_f = casadi.Function('sos_f', {Qvar}, {Hf Df Qobj}, struct('allow_free', true));
sosprob_g = casadi.Function('sos_g', {Qvar}, {Dg Qcon}, struct('allow_free', true));

% map sum-of-squares decision variables to matrix variables
map_x = blkdiag(speye(nnz_lin_x), Mp_x);
Qvar_mapped = map_x*Qvar_sdp;

% evaluate functions
[sdp_Hf, sdp_Jf, sdp_f] = sosprob_f(Qvar_mapped);
[sdp_Jg, sdp_g] = sosprob_g(Qvar_mapped);

% map sum-of-squares constraints to matrix variables
map_g = blkdiag(sparse(nnz_lin_g,0), Mp_g);

% build SDP problem
sdp.x = permMat*[Qvar_sdp; Qcon_G];
sdp.f = sdp_f;
sdp.g = sdp_g - map_g*Qcon_G;
sdp.p = Qpar;

% store derivatives in SDP problem
permMat_t = permMat';
n_vars = nnz_lin_x + nnz_gram_x;

% store derivatives
sdp.derivatives.Hf = permMat_t(1:n_vars, :)'*map_x'*sdp_Hf*map_x*permMat_t(1:(n_vars),:);
sdp.derivatives.Jf = sdp_Jf*map_x*permMat_t(1:n_vars, :);
sdp.derivatives.Jg = sdp_Jg*map_x*permMat_t(1:n_vars, :)-map_g*permMat_t(n_vars+1:end, :);

%sdp.derivatives.Hf = blockcat(map_x'*sdp_Hf*map_x, ...
%                              sparse(nnz_lin_x+nnz_gram_x,nnz_gram_g), ...
%                              sparse(nnz_gram_g,nnz_lin_x+nnz_gram_x), ...
%                              sparse(nnz_gram_g,nnz_gram_g));
%sdp.derivatives.Jf = horzcat(sdp_Jf*map_x, sparse(1,nnz_gram_g));
%sdp.derivatives.Jg = horzcat(sdp_Jg*map_x, -map_g);

% SDP options
sdpopt = opts.sdpsol_options;

% define the cones in the sdp level
sdpopt.Kx.lin = nnz_lin_x;

% split the matrices using mat2cell
cell_Kmat_x = mat2cell(Ksdp_x_s, [Ns Nds Nsds], 1);
cell_Kmat_g = mat2cell(Ksdp_g_s, [Ms Mds Msds], 1);

% split and assign to struct fields
% PSD cones
if ~isempty(cell_Kmat_x{1}) || ~isempty(cell_Kmat_g{1})
    sdpopt.Kx.psd = [cell_Kmat_x{1}; cell_Kmat_g{1}];
end
% DD cones
if ~isempty(cell_Kmat_x{2}) || ~isempty(cell_Kmat_g{2})
    sdpopt.Kx.dd  = [cell_Kmat_x{2}; cell_Kmat_g{2}];
end
% SDD cones
if ~isempty(cell_Kmat_x{3}) || ~isempty(cell_Kmat_g{3})
    sdpopt.Kx.sdd = [cell_Kmat_x{3}; cell_Kmat_g{3}];
end

% assign number of linear cones in constraints
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

% % coordinates of SOS solution
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

