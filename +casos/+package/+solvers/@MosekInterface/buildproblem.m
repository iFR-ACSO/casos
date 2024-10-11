function buildproblem(obj)
% Construct MOSEK problem description from problem structure.
%
% MOSEK solves semidefinite problems in the form
%
%   min tr(Cbar*X) + c'*x
%   st. lba <= tr(Abar*X) + A*x <= uba,
%       tr(Fbar*X) + F*x + g in Kg,
%       lbx <= x <= ubx, x in Kx,
%       X is psd
%
% where Kx is a variable cone (linear, quadratic, and more) and Kg is an
% affine constraint cone (including vectorized LMIs).
%
% The generic conic problem in casos has the form
%
%   min 1/2*x'*h*x + g'*x
%   st. a*x in Kc, x in Kx,
%
% with Lagrange multipliers satisfying
%
%   h*x + g + a'*lam_a + lam_x = 0
%   lam_a in -Kc*, lam_x in -Kx*
%
% where Kc and Kx are cones composed of
%   - box constraints [lb ub] (l)
%   - Lorentz (quadratic) cone (q)
%   - rotated Lorentz cone (r)
%   - cone of PSD matrices (s)
%
% and Kc* and Kx* are the dual cones of Kc and Kx, respectively.
%
% Lorentz cone, rotated Lorentz cone, and PSD cone can be shifted by a
% lower bound (cb).

opts = obj.opts;

% retrieve MOSEK symbolic constants
[rcode,res] = mosekopt('symbcon echo(0)');

assert(rcode == 0,'Call to mosekopt failed: %s (%s).',res.rmsg,res.rcodestr);
symbcon = res.symbcon;

% symbolic variables
h = obj.args_in.h;
g = obj.args_in.g;
a = obj.args_in.a;

% symbolic bounds
lba = obj.args_in.lba;
uba = obj.args_in.uba;
cba = obj.args_in.cba;
lbx = obj.args_in.lbx;
ubx = obj.args_in.ubx;
cbx = obj.args_in.cbx;

% problem size
[m,n] = size(a);

% obtain cones
Kx = opts.Kx;
Kc = opts.Kc;

% number of variables per cone type
Nx.lin = (obj.getdimc(Kx,'lin'));
Nx.lor = (obj.getdimc(Kx,'lor'));
Nx.rot = (obj.getdimc(Kx,'rot'));
Nx.psd = (obj.getdimc(Kx,'psd'));
Nx.pow = (obj.getdimc(Kx,'pow'));
Nx.dpow = (obj.getdimc(Kx,'dpow'));
Nx.exp = (obj.getdimc(Kx,'exp'));
Nx.dexp = (obj.getdimc(Kx,'dexp'));

% number of constraints per cone type
Na.lin = (obj.getdimc(Kc,'lin'));
Na.lor = (obj.getdimc(Kc,'lor'));
Na.rot = (obj.getdimc(Kc,'rot'));
Na.psd = (obj.getdimc(Kc,'psd'));
Na.pow = (obj.getdimc(Kc,'pow'));
Na.dpow = (obj.getdimc(Kc,'dpow'));
Na.exp = (obj.getdimc(Kc,'exp'));
Na.dexp = (obj.getdimc(Kc,'dexp'));

% number of variables in the analytic cones
Nx_p = Nx.exp + Nx.dexp + obj.getnumc(Kx,'pow') + obj.getnumc(Kx,'dpow');
% number of vector-valued variables
Nx_v = n - sum(Nx.psd.^2) - Nx_p;
% number of quadratic variables
Nx_q = sum([Nx.lor Nx.rot]);
% number of conic constraints
Na_c = m - Na.lin;
% number of quadratic constraints
Na_q = sum([Na.lor Na.rot]);
% number of analytic constraints
Na_p = Na.exp + Na.dexp + obj.getnumc(Kc,'pow') + obj.getnumc(Kc,'dpow');

% number of variable constraints
% to write into affine cone constraints
Nx_acc = Nx_q + Nx_p;

% semidefinite variables (vectorized)
Nx_S = Nx.psd.*(Nx.psd+1)/2;
Na_S = Na.psd.*(Na.psd+1)/2;
% affine cone constraints (vectorized)
Na_C = Na_c + sum(Na_S - Na.psd.^2);

% separate conic bound for vector, matrix, and analytic conic variables
[cbx_v,cbx_s,cbx_p] = separate(cbx,[Nx_q sum(Nx.psd.^2) Nx_p]);

% separate linear cost for vector, matrix, and analytic variables
% C' = | c : Cbar |
[Clin,Cbar,Cpoly] = separate(g,[Nx_v sum(Nx.psd.^2) Nx_p],1);
% linear cost vector (vector and analytic variables)
prob.c = [Clin; Cpoly];
% symmetric cost matrices Cbar_j as stacked vectorization
% Cbar = [Cbar1(:); ...; CbarN(:)]
% get nonzero elements and subindices (j,k,l)
[barv.c,barc.subj,~,barc.subk,barc.subl] = obj.sdp_vec(Cbar,Nx.psd,1);

% separate linear constraints and affine cone constraints
%     | a : Abar | 
% A = |----------|
%     |    F     |
[A,F] = separate(a,[Na.lin Na_c],n);

% separate linear constraints for vector, matrix, and analytic variables
[Alin,Abar,Apoly] = separate(A,Na.lin,[Nx_v sum(Nx.psd.^2) Nx_p]);
% linear constraints matrix (vector and analytic variables)
prob.a = [Alin Apoly];
% symmetric constraint matrices Abar_ij as stacked vectorizations
%        | Abar11(:)' ... Abar1N(:)' |
% Abar = |     :       :      :      |
%        | AbarM1(:)' ... AbarMN(:)' |
% get nonzero elements and subindices (i,j,k,l)
[barv.a,bara.subi,bara.subj,bara.subk,bara.subl] = obj.sdp_vec(Abar,Nx.psd,1,2);
% rewrite 
%   lba <= Abar*X <= uba, X in Kx + cbx
% into
%   lba <= Abar*(X'+cbx) <= uba, X' in Kx
Acb = Abar*cbx_s;
% linear bounds
prob.blc = lba - Acb;
prob.buc = uba - Acb;

% vectorize semidefinite domain for affine cone constraints
%     | Fvec |
%     |------|
% F = | Fmat |
%     |------|
%     | Fpol |
Fc = mat2cell(   F,[Na_q sum(Na.psd.^2) Na_p],n);
gc = mat2cell(-cba,[Na_q sum(Na.psd.^2) Na_p],1);
% vector-based vectorization
Fmat = obj.sdp_vec(Fc{2},Na.psd,[],1);
gmat = obj.sdp_vec(gc{2},Na.psd,[],1);
% build affine cone constraints
Facc = vertcat(Fc{1},Fmat,Fc{3});
gacc = vertcat(gc{1},gmat,gc{3});

% separate affine cone constraints for vector, matrix, and analytic variables
% F = | f : Fbar |
[Flin,Fbar,Fpoly] = separate(Facc,Na_C,[Nx_v sum(Nx.psd.^2) Nx_p]);
% affine cone constraint matrix
prob.f = [
    Flin Fpoly
    casadi.DM.triplet((1:Nx_acc)-1,Nx.lin+(1:Nx_acc)-1,ones(Nx_acc,1),Nx_acc,Nx_v+Nx_p)
];
% symmetric constraint matrices as stacked vectorization
% get nonzero elements and subindices (i,j,k,l)
[barv.f,barf.subi,barf.subj,barf.subk,barf.subl] = obj.sdp_vec(Fbar,Nx.psd,1,2);
% rewrite 
%   Fbar*X + g in Kc, X in Kx + cbx
% into
%   Fbar*(X'+cbx) + g in Kc, X' in Kx
Fcb = Fbar*cbx_s;
% constant cone
prob.g = [gacc + Fcb; cbx_v; cbx_p];

% linear and vector state constraints
prob.blx = [lbx; -inf(Nx_acc,1)];
prob.bux = [ubx; +inf(Nx_acc,1)];

% handle quadratic cost function
% MOSEK cannot solve conic problems with quadratic cost
if nnz(h) > 0
    % only SX supports Cholesky decomposition
    H = casadi.SX.sym('H',sparsity(h));
    chol_f = casadi.Function('chol',{H},{chol(H)});
    % rewrite quadratic cost
    %
    %   min_x 1/2 x'*Q*x + c'*x
    %
    % into 
    %
    %   min_{x,y} c'*x + y, s.t. 1 + y >= ||(sqrt(2)*U*x, 1 - y)||
    %
    % with additional variable y and Cholesky decomposition U'*U = Q
    U = chol_f(h);
    % build affine cone constraint 
    % L(y,x) + k = (1+y, sqrt(2)*U*x, 1-y) in SOC
    % note: additional variable y is first decision variable
    L = [1 casadi.DM(1,n); casadi.DM(n,1) sqrt(2)*U; -1 casadi.DM(1,n)];
    k = [1; casadi.DM(n,1); 1];
    % number of additional variables and constraints
    Nx_cost = 1;
    Na_cost = n + 2;
    % separate ACC for vector, matrix, and analytic variables
    [Llin,Lbar,Lpoly] = separate(L,Na_cost,[Nx_v+1 sum(Nx.psd.^2) Nx_p]);
    % get nonzero elements and subindices (i,j,k,l) for Lbar
    [barv_l,barl.subi,barl.subj,barl.subk,barl.subl] = obj.sdp_vec(Lbar,Nx.psd,1,2);
    % add to regular (linear) cost
    prob.c = [1; prob.c];
    % add to regular linear constraints
    prob.a = [casadi.DM(Na.lin,1) prob.a];
    % add to regular affine cone constraints
    prob.f = [casadi.DM(Na_C+Nx_acc,1) prob.f; Llin Lpoly];
    prob.g = [prob.g; k];
    barv.f = [barv.f barv_l];
    barf.subi = [barf.subi (Na_C+Nx_acc)+barl.subi];
    barf.subj = [barf.subj barl.subj];
    barf.subk = [barf.subk barl.subk];
    barf.subl = [barf.subl barl.subl];
    % add to regular decision variables (y >= 0)
    prob.blx = [  0; prob.blx];
    prob.bux = [inf; prob.bux];
    % add to conic domain
    acc_cost = [symbcon.MSK_DOMAIN_QUADRATIC_CONE Na_cost];

else
    Nx_cost = 0;
    Na_cost = 0;
    acc_cost = [];
end

% options for Casadi functions
% NOTE: As per Casadi v3.6.5, functions' input and output names 
% must be mutually exclusive unless explicity permitted.
fopt = struct('allow_duplicate_io_names',true);
% return MOSEK prob structure
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),struct2cell(prob),fieldnames(obj.args_in),fieldnames(prob),fopt);
% return bar values
obj.barv = casadi.Function('v',struct2cell(obj.args_in),struct2cell(barv),fieldnames(obj.args_in),fieldnames(barv),fopt);

% build conic information
Accs = [
    arrayfun(@(l) [symbcon.MSK_DOMAIN_QUADRATIC_CONE  l], Na.lor(:), 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_RQUADRATIC_CONE l], Na.rot(:), 'UniformOutput',false)
    arrayfun(@(d) [symbcon.MSK_DOMAIN_SVEC_PSD_CONE   d], Na_S(:), 'UniformOutput',false)
    arrayfun(@(s) [symbcon.MSK_DOMAIN_PRIMAL_POWER_CONE s.dim length(s.alpha) s.alpha(:)'], Na.pow(:), 'UniformOutput',false)
    arrayfun(@(s) [symbcon.MSK_DOMAIN_DUAL_POWER_CONE s.dim length(s.alpha) s.alpha(:)'], Na.dpow(:), 'UniformOutput',false)
          repmat({[symbcon.MSK_DOMAIN_PRIMAL_EXP_CONE 3]}, Na.exp/3, 1)
          repmat({[symbcon.MSK_DOMAIN_DUAL_EXP_CONE   3]}, Na.dexp/3, 1)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_QUADRATIC_CONE  l], Nx.lor(:), 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_RQUADRATIC_CONE l], Nx.rot(:), 'UniformOutput',false)
    arrayfun(@(s) [symbcon.MSK_DOMAIN_PRIMAL_POWER_CONE s.dim length(s.alpha) s.alpha(:)'], Nx.pow(:), 'UniformOutput',false)
    arrayfun(@(s) [symbcon.MSK_DOMAIN_DUAL_POWER_CONE s.dim length(s.alpha) s.alpha(:)'], Nx.dpow(:), 'UniformOutput',false)
          repmat({[symbcon.MSK_DOMAIN_PRIMAL_EXP_CONE 3]}, Nx.exp/3, 1)
          repmat({[symbcon.MSK_DOMAIN_DUAL_EXP_CONE   3]}, Nx.dexp/3, 1)
    acc_cost
];
cone.bardim = Nx.psd;
cone.barc = barc;
cone.bara = bara;
cone.barf = barf;
cone.accs = horzcat(Accs{:});
% store static information
obj.cone = cone;

% parse MOSEK solution into (x,cost,lam_a,lam_x)
sol.pobjval = casadi.MX.sym('pobjval');
sol.xx   = casadi.MX.sym('xx',[Nx_cost+Nx_v+Nx_p 1]);
sol.barx = casadi.MX.sym('barx',[sum(Nx_S) 1]);
sol.slc  = casadi.MX.sym('slc',[Na.lin 1]);
sol.suc  = casadi.MX.sym('suc',[Na.lin 1]);
sol.slx  = casadi.MX.sym('slx',[Nx_cost+Nx.lin+Nx_q+Nx_p 1]);
sol.sux  = casadi.MX.sym('sux',[Nx_cost+Nx.lin+Nx_q+Nx_p 1]);
sol.doty = casadi.MX.sym('doty',[Na_q+sum(Na_S)+Na_p+Nx_q+Nx_p+Na_cost 1]);
sol.bars = casadi.MX.sym('bars',[sum(Nx_S) 1]);

% dual variables corresponding to linear variables
[~,Slx,~] = separate(sol.slx,[Nx_cost Nx.lin Nx_q Nx_p],1);
[~,Slu,~] = separate(sol.sux,[Nx_cost Nx.lin Nx_q Nx_p],1);
% dual variables corresponding to affine conic constraints
[Yaq,Yas,Yap,Yxq,Yxp,Ycost] = separate(sol.doty,[Na_q sum(Na_S) Na_p Nx_q Nx_p Na_cost],1);
% de-vectorize SDP primal and dual variables (no scaling)
Xc_s = obj.sdp_mat(sol.barx,Nx.psd,1) + cbx_s;
Sc_s = obj.sdp_mat(sol.bars,Nx.psd,1);
% de-vectorize duals corresponding to semidefinite constraints
Yc_s = obj.sdp_mat(Yas,Na.psd,1);
% multipliers for box constraints
lam_a_l = sol.suc - sol.slc;
lam_x_l = Slu - Slx;
% multipliers for quadratic constraints
lam_a_q = -Yaq;
lam_x_q = -Yxq;
% multipliers for SDP constraints
lam_a_s = -vertcat(Yc_s);
lam_x_s = -vertcat(Sc_s);
% multipliers for analytic constraints
lam_a_p = -Yap;
lam_x_p = -Yxp;
% build multipliers
lam_a = [lam_a_l; lam_a_q; lam_a_s; lam_a_p];
lam_x = [lam_x_l; lam_x_q; lam_x_s; lam_x_p];
% build solution
[~,Xc_l,Xc_p] = separate(sol.xx,[Nx_cost Nx_v Nx_p],1);
sol_x = vertcat(Xc_l,Xc_s,Xc_p);
% cost
cost = sol.pobjval;

obj.ghan = casadi.Function('g',[struct2cell(sol); struct2cell(obj.args_in)],{sol_x cost lam_a lam_x},[fieldnames(sol); fieldnames(obj.args_in)],obj.names_out);

end

function varargout = separate(A,varargin)
% Separate array into subarrays.

    varargout = mat2cell(A,varargin{:});
end
