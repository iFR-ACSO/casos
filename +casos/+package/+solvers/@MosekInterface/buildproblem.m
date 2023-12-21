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
Nx.l = (obj.getdimc(Kx,'l'));
Nx.q = (obj.getdimc(Kx,'q'));
Nx.r = (obj.getdimc(Kx,'r'));
Nx.s = (obj.getdimc(Kx,'s'));

% number of constraints per cone type
Na.l = (obj.getdimc(Kc,'l'));
Na.q = (obj.getdimc(Kc,'q'));
Na.r = (obj.getdimc(Kc,'r'));
Na.s = (obj.getdimc(Kc,'s'));

% assert(Nx.q+Nx.r == 0, 'Quadratic variables not supporte yet.')
% assert(Na.q+Na.q == 0, 'Quadratic constraints not supporte yet.')

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

% separate linear cost for vector and matrix variables
Cc = mat2cell(g,[Nx_v sum(Nx.s.^2)],1);
% linear cost vector
prob.c = Cc{1};
% symmetric cost matrices Cbar_j as stacked vectorization
% Cbar = [Cbar1(:); ...; CbarN(:)]
Cbar = Cc{2};
% get nonzero elements and subindices (j,k,l)
[barv.c,barc.subj,~,barc.subk,barc.subl] = obj.sdp_vec(Cbar,Nx.s,1);
% % get nonzero indices w.r.t. Cbar
% Ic = find(sparsity(Cbar));
% % get subindices (j,k,l) and linear indices for nonzero elements (tril)
% [barc.subj,barc.subk,barc.subl,Ic] = get_barsub(Nx.s,Ic);
% % store nonzeros elements
% barv.c = Cbar(Ic); %vertcat(barc_val{Ic});

% separate linear constraints and affine cone constraints
%     | a : Abar | 
% A = |----------|
%     |    F     |
Ac = mat2cell(a,[Na.l Na_c],n);

% separate linear constraints for vector and matrix variables
Alin = mat2cell(Ac{1},Na.l,[Nx_v sum(Nx.s.^2)]);
% linear constraints matrix
prob.a = Alin{1};
% symmetric constraint matrices Abar_ij as stacked vectorizations
%        | Abar11(:)' ... Abar1N(:)' |
% Abar = |     :       :      :      |
%        | AbarM1(:)' ... AbarMN(:)' |
Abar = Alin{2};
% get nonzero elements and subindices (i,j,k,l)
[barv.a,bara.subi,bara.subj,bara.subk,bara.subl] = obj.sdp_vec(Abar,Nx.s,1,2);
% % get nonzero indices w.r.t. Abar
% [iA,jA] = ind2sub(size(Abar),find(sparsity(Abar)));
% % get subindices (j,k,l) and linear indices for nonzero elements (tril)
% [bara.subj,bara.subk,bara.subl,jA,iA] = get_barsub(Nx.s,jA,iA);
% % subindices i for nonzero elements
% bara.subi = iA;
% % store nonzeros elements
% Ia = sub2ind(size(Abar),iA,jA);
% barv.a = Abar(Ia); %vertcat(bara_val{Ia,Ja});
% linear bounds
prob.blc = lba;
prob.buc = uba;

% vectorize semidefinite domain for affine cone constraints
%     | Fvec |
% F = |------|
%     | Fmat |
Fc = mat2cell(Ac{2},[Na_q sum(Na.s.^2)],n);
gc = mat2cell( -cba,[Na_q sum(Na.s.^2)],1);
% vector-based vectorization
Fmat = obj.sdp_vec(Fc{2},Na.s,[],1);
gmat = obj.sdp_vec(gc{2},Na.s,[],1);
% Fmat = cellfun(@obj.sdp_vec, mat2cell(Fc{2},Na.s.^2,n), 'UniformOutput',false);
% gmat = cellfun(@obj.sdp_vec, mat2cell(gc{2},Na.s.^2,1), 'UniformOutput',false);
% build affine cone constraints
Facc = vertcat(Fc{1},Fmat);
gacc = vertcat(gc{1},gmat);

% separate affine cone constraints for vector and matrix variables
% F = | f : Fbar |
Fcc = mat2cell(Facc,Na_c,[Nx_v sum(Nx.s.^2)]);
% affine cone constraint matrix
prob.f = [
    Fcc{1}
    sparse(1:Nx_q,Nx.l+(1:Nx_q),ones(Nx_q,1),Nx_q,Nx_v)
];
% constant cone
prob.g = [gacc; zeros(Nx_q,1)];
% symmetric constraint matrices as stacked vectorization
Fbar = Fcc{2};
% get nonzero elements and subindices (i,j,k,l)
[barv.f,barf.subi,barf.subj,barf.subk,barf.subl] = obj.sdp_vec(Fbar,Nx.s,1,2);
% % get nonzero indices w.r.t. Fbar
% [iF,jF] = ind2sub(size(Fbar),find(sparsity(Fbar)));
% % get subindices (j,k,l) and linear indices for nonzero elements (tril)
% [barf.subj,barf.subk,barf.subl,jF,iF] = get_barsub(Nx.s,jF,iF);
% % subindices i for nonzero elements
% barf.subi = iF;
% % store nonzeros elements
% If = sub2ind(size(Fbar),iF,jF);
% barv.f = Fbar(If); %vertcat(barf_val{If,Jf});

% linear and vector state constraints
prob.blx = [lbx; -inf(Nx_q,1)];
prob.bux = [ubx; +inf(Nx_q,1)];

% return MOSEK prob structure
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),struct2cell(prob),fieldnames(obj.args_in),fieldnames(prob));
% return bar values
obj.barv = casadi.Function('v',struct2cell(obj.args_in),struct2cell(barv),fieldnames(obj.args_in),fieldnames(barv));

% build conic information
Accs = [
    arrayfun(@(l) [symbcon.MSK_DOMAIN_QUADRATIC_CONE  l], Na.q, 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_RQUADRATIC_CONE l], Na.r, 'UniformOutput',false)
    arrayfun(@(d) [symbcon.MSK_DOMAIN_SVEC_SDP_CONE   d], Na.s, 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_QUADRATIC_CONE  l], Nx.q, 'UniformOutput',false)
    arrayfun(@(l) [symbcon.MSK_DOMAIN_RQUADRATIC_CONE l], Nx.r, 'UniformOutput',false)
];
cone.bardim = Kx.s;
cone.barc = barc;
cone.bara = bara;
cone.barf = barf;
cone.accs = horzcat(Accs{:});
% store static information
obj.cone = cone;

% parse MOSEK solution into (x,cost,lam_a,lam_x)
sol.pobjval = casadi.MX.sym('pobjval');
sol.xx   = casadi.MX.sym('xx',[Nx_v 1]);
sol.barx = casadi.MX.sym('barx',[sum(Nx_S) 1]);
% sol.barx = arrayfun(@(d) casadi.MX.sym('Xbar',[d 1]), Nx_S, 'UniformOutput',false);
sol.slc  = casadi.MX.sym('slc',[Na.l 1]);
sol.suc  = casadi.MX.sym('suc',[Na.l 1]);
sol.slx  = casadi.MX.sym('slx',[Nx.l+Nx_q 1]);
sol.sux  = casadi.MX.sym('sux',[Nx.l+Nx_q 1]);
sol.doty = casadi.MX.sym('doty',[Na_q+sum(Na_S)+Nx_q 1]);
sol.bars = casadi.MX.sym('bars',[sum(Nx_S) 1]);
% sol.bars = arrayfun(@(d) casadi.MX.sym('Sbar',[d 1]), Nx_S, 'UniformOutput',false);

% dual variables corresponding to affine conic constraints
Yc = mat2cell(sol.doty,[Na_q sum(Na_S) Nx_q],1);
% de-vectorize SDP primal and dual variables (no scaling)
Xc_s = obj.sdp_mat(sol.barx,Nx.s,1);
Sc_s = obj.sdp_mat(sol.bars,Nx.s,1);
% Xc_s = cellfun(@(X) obj.sdp_mat(X,1), mat2cell(sol.barx,Nx_S,1), 'UniformOutput', false);
% Sc_s = cellfun(@(S) obj.sdp_mat(S,1), mat2cell(sol.bars,Nx_S,1), 'UniformOutput', false);
% de-vectorize duals corresponding to semidefinite constraints
Yc_s = obj.sdp_mat(Yc{2},Na.s,1);
% Yc_s = cellfun(@obj.sdp_mat, mat2cell(Yc{2},Na_S,1), 'UniformOutput',false);
% multipliers for box constraints
lam_a_l = sol.suc - sol.slc;
lam_x_l = sol.sux(1:Nx.l) - sol.slx(1:Nx.l);
% multipliers for quadratic constraints
lam_a_q = -Yc{1};
lam_x_q = -Yc{3};
% multipliers for SDP constraints
lam_a_s = -vertcat(Yc_s);
lam_x_s = -vertcat(Sc_s);
% build multipliers
lam_a = [lam_a_l; lam_a_q; lam_a_s];
lam_x = [lam_x_l; lam_x_q; lam_x_s];
% build solution
sol_x = vertcat(sol.xx,Xc_s);
% cost
cost = sol.pobjval;

obj.ghan = casadi.Function('g',struct2cell(sol),{sol_x cost lam_a lam_x},fieldnames(sol),obj.names_out);

end

function [j,k,l,I,J] = get_barsub(s,I,J)
% Return subindices for into nonzeros elements of stacked coefficients.

    I = reshape(I,1,[]); % ensure row vector of linear indices.
    s = reshape(s,1,[]); % ensure row vector of matrix dimensions.

    % number of matrix components before j-th matrix variables
    S = cumsum([0 s(1:end-1).^2]);

    % compute index j of corresponding matrix variable
    j = sum(I > S', 1); % MOSEK interface is 1-based

    % compute linear indices of elements in each matrix
    I0 = I - S(j);

    % compute indices (k,l) of elements in matrix
    % as (remainder,quotient) of linear indices divided by matrix size
    l = ceil(I0 ./ s(j)); % col
    k = I0 - s(j).*(l-1); % row

    % only keep lower-triangular elements
    Itril = (k >= l);
    % remove indices for upper triangular
    I(~Itril) = [];
    j(~Itril) = [];
    l(~Itril) = [];
    k(~Itril) = [];

    % optional row indices
    if nargin > 2
        J(~Itril) = [];
    end
end
