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

if ~isfield(obj.opts,'cholesky_method'), obj.opts.cholesky_method = 'numerical'; end

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
Nx.l = (obj.getdimc(Kx,'lin'));
Nx.q = (obj.getdimc(Kx,'lor'));
Nx.r = (obj.getdimc(Kx,'rot'));
Nx.s = (obj.getdimc(Kx,'psd'));

% number of constraints per cone type
Na.l = (obj.getdimc(Kc,'lin'));
Na.q = (obj.getdimc(Kc,'lor'));
Na.r = (obj.getdimc(Kc,'rot'));
Na.s = (obj.getdimc(Kc,'psd'));

% number of vector-valued variables
Nx_v = n - sum(Nx.s.^2);
% number of quadratic variables
Nx_q = sum([Nx.q Nx.r]);
% number of conic constraints
Na_c = m - Na.l;
% number of quadratic constraints
Na_q = sum([Na.q Na.r]);

% semidefinite variables (vectorized)
Nx_S = Nx.s.*(Nx.s+1)/2;
Na_S = Na.s.*(Na.s+1)/2;
% affine cone constraints (vectorized)
Na_C = Na_c + sum(Na_S - Na.s.^2);

% separate conic bound for vector and matrix conic variables
[cbx_v,cbx_s] = separate(cbx,[Nx_q sum(Nx.s.^2)]);

% separate linear cost for vector and matrix variables
% C' = | c : Cbar |
[Clin,Cbar] = separate(g,[Nx_v sum(Nx.s.^2)],1);
% linear cost vector
prob.c = Clin;
% symmetric cost matrices Cbar_j as stacked vectorization
% Cbar = [Cbar1(:); ...; CbarN(:)]
% get nonzero elements and subindices (j,k,l)
[barv.c,barc.subj,~,barc.subk,barc.subl] = obj.sdp_vec(Cbar,Nx.s,1);

% separate linear constraints and affine cone constraints
%     | a : Abar | 
% A = |----------|
%     |    F     |
[A,F] = separate(a,[Na.l Na_c],n);

% separate linear constraints for vector and matrix variables
[Alin,Abar] = separate(A,Na.l,[Nx_v sum(Nx.s.^2)]);
% linear constraints matrix
prob.a = Alin;
% symmetric constraint matrices Abar_ij as stacked vectorizations
%        | Abar11(:)' ... Abar1N(:)' |
% Abar = |     :       :      :      |
%        | AbarM1(:)' ... AbarMN(:)' |
% get nonzero elements and subindices (i,j,k,l)
[barv.a,bara.subi,bara.subj,bara.subk,bara.subl] = obj.sdp_vec(Abar,Nx.s,1,2);
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
% F = |------|
%     | Fmat |
Fc = mat2cell(   F,[Na_q sum(Na.s.^2)],n);
gc = mat2cell(-cba,[Na_q sum(Na.s.^2)],1);
% vector-based vectorization
Fmat = obj.sdp_vec(Fc{2},Na.s,[],1);
gmat = obj.sdp_vec(gc{2},Na.s,[],1);
% build affine cone constraints
Facc = vertcat(Fc{1},Fmat);
gacc = vertcat(gc{1},gmat);

% separate affine cone constraints for vector and matrix variables
% F = | f : Fbar |
[Flin,Fbar] = separate(Facc,Na_C,[Nx_v sum(Nx.s.^2)]);
% affine cone constraint matrix
prob.f = [
    Flin
    casadi.DM.triplet((1:Nx_q)-1,Nx.l+(1:Nx_q)-1,ones(Nx_q,1),Nx_q,Nx_v)
];
% symmetric constraint matrices as stacked vectorization
% get nonzero elements and subindices (i,j,k,l)
[barv.f,barf.subi,barf.subj,barf.subk,barf.subl] = obj.sdp_vec(Fbar,Nx.s,1,2);
% rewrite 
%   Fbar*X + g in Kc, X in Kx + cbx
% into
%   Fbar*(X'+cbx) + g in Kc, X' in Kx
Fcb = Fbar*cbx_s;
% constant cone
prob.g = [gacc + Fcb; cbx_v];

% linear and vector state constraints
prob.blx = [lbx; -inf(Nx_q,1)];
prob.bux = [ubx; +inf(Nx_q,1)];


% arguments to problem
args_in = obj.args_in;

% handle quadratic cost function
% MOSEK cannot solve conic problems with quadratic cost

if nnz(h) > 0
    
    % use casadi to compute cholesky
    if strcmp(obj.opts.cholesky_method,'analytical') 

        % only SX supports Cholesky decomposition
        H = casadi.SX.sym('H',sparsity(h));
        %Performs an LDL transformation [L,D] = ldl(A) and returns 
        % diag(sqrt(D))*L'
        chol_f = casadi.Function('chol',{H},{chol(H)});
    
        % ldl_f = casadi.Function('ldl',{H},{ldl(H)});
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
    
    else % parameterize solver with cholesky sparsity pattern and compute numerically online

        % get sparisty pattern from h
        H = sparsity(h);
        [rr,~] = get_triplet(H);
        
        i = min(rr); j = max(rr); % row numbers are sorted
        
        spU = casadi.Sparsity.upper(j - i + 1);
        
        spU.enlarge(size(H,1), size(H,2), i:j, i:j);
        
        U = casadi.MX.sym('U', spU); %+casadi.Sparsity.diag(length(H)));

       % evaluate Cholesky decomposition in situ
       args_in.h = U;
   
    end

    % build affine cone constraint 
    % L(y,x) + k = (1+y, sqrt(2)*U*x, 1-y) in SOC
    % note: additional variable y is first decision variable
    L = [1 casadi.DM(1,n); casadi.DM(n,1) sqrt(2)*U; -1 casadi.DM(1,n)];
    k = [1; casadi.DM(n,1); 1];


    % number of additional variables and constraints
    Nx_cost = 1;
    Na_cost = n + 2;
    % separate ACC for vector and matrix variables
    [Llin,Lbar] = separate(L,Na_cost,[Nx_v+1 sum(Nx.s.^2)]);
    % get nonzero elements and subindices (i,j,k,l) for Lbar
    [barv_l,barl.subi,barl.subj,barl.subk,barl.subl] = obj.sdp_vec(Lbar,Nx.s,1,2);
    % add to regular (linear) cost
    prob.c = [1; prob.c];
    % add to regular linear constraints
    prob.a = [casadi.DM(Na.l,1) prob.a];
    % add to regular affine cone constraints
    prob.f = [casadi.DM(Na_C+Nx_q,1) prob.f; Llin];
    prob.g = [prob.g; k];
    barv.f = [barv.f barv_l];
    barf.subi = [barf.subi (Na_C+Nx_q)+barl.subi];
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
obj.fhan = casadi.Function('f',struct2cell(args_in),struct2cell(prob),fieldnames(args_in),fieldnames(prob),fopt);
% return bar values
obj.barv = casadi.Function('v',struct2cell(args_in),struct2cell(barv),fieldnames(args_in),fieldnames(barv),fopt);

% build conic information
Accs = [
    arrayfun(@(l) [symbcon.MSK_DOMAIN_QUADRATIC_CONE  l], Na.q, 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_RQUADRATIC_CONE l], Na.r, 'UniformOutput',false)
    arrayfun(@(d) [symbcon.MSK_DOMAIN_SVEC_PSD_CONE   d], Na_S, 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_QUADRATIC_CONE  l], Nx.q, 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_RQUADRATIC_CONE l], Nx.r, 'UniformOutput',false)
    acc_cost
];
cone.bardim = Nx.s;
cone.barc = barc;
cone.bara = bara;
cone.barf = barf;
cone.accs = horzcat(Accs{:});
% store static information
obj.cone = cone;

% parse MOSEK solution into (x,cost,lam_a,lam_x)
sol.pobjval = casadi.MX.sym('pobjval');
sol.xx   = casadi.MX.sym('xx',[Nx_cost+Nx_v 1]);
sol.barx = casadi.MX.sym('barx',[sum(Nx_S) 1]);
sol.slc  = casadi.MX.sym('slc',[Na.l 1]);
sol.suc  = casadi.MX.sym('suc',[Na.l 1]);
sol.slx  = casadi.MX.sym('slx',[Nx_cost+Nx.l+Nx_q 1]);
sol.sux  = casadi.MX.sym('sux',[Nx_cost+Nx.l+Nx_q 1]);
sol.doty = casadi.MX.sym('doty',[Na_q+sum(Na_S)+Nx_q+Na_cost 1]);
sol.bars = casadi.MX.sym('bars',[sum(Nx_S) 1]);

% dual variables corresponding to linear variables
[~,Slx,~] = separate(sol.slx,[Nx_cost Nx.l Nx_q],1);
[~,Slu,~] = separate(sol.sux,[Nx_cost Nx.l Nx_q],1);
% dual variables corresponding to affine conic constraints
[Yaq,Yas,Yxq,Ycost] = separate(sol.doty,[Na_q sum(Na_S) Nx_q Na_cost],1);
% de-vectorize SDP primal and dual variables (no scaling)
Xc_s = obj.sdp_mat(sol.barx,Nx.s,1) + cbx_s;
Sc_s = obj.sdp_mat(sol.bars,Nx.s,1);
% de-vectorize duals corresponding to semidefinite constraints
Yc_s = obj.sdp_mat(Yas,Na.s,[]);
% multipliers for box constraints
lam_a_l = sol.suc - sol.slc;
lam_x_l = Slu - Slx;
% multipliers for quadratic constraints
lam_a_q = -Yaq;
lam_x_q = -Yxq;
% multipliers for SDP constraints
lam_a_s = -vertcat(Yc_s);
lam_x_s = -vertcat(Sc_s);
% build multipliers
lam_a = [lam_a_l; lam_a_q; lam_a_s];
lam_x = [lam_x_l; lam_x_q; lam_x_s];
% build solution
[~,Xc_l] = separate(sol.xx,[Nx_cost Nx_v],1);
sol_x = vertcat(Xc_l,Xc_s);
% cost
cost = sol.pobjval;

obj.ghan = casadi.Function('g',[struct2cell(sol); struct2cell(obj.args_in)],{sol_x cost lam_a lam_x},[fieldnames(sol); fieldnames(obj.args_in)],obj.names_out);


% fill info struct
obj.solver_info.name        = 'mosek';
obj.solver_info.n_decVar    = length(obj.args_in.x0);
obj.solver_info.size_A      = size(obj.args_in.a);

end

function varargout = separate(A,varargin)
% Separate array into subarrays.

    varargout = mat2cell(A,varargin{:});
end
