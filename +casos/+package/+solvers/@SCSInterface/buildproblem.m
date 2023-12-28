function buildproblem(obj)
% Construct SCS problem description (P,A,b,c,K) from problem structure.
%
% SCS solves semidefinite problems in the form
%
%   min 1/2*x'*P*x + c'*x 
%   s.t. A*x + s = b, s in K
%
% with dual variabes satisfying
%
%   P*x + A'*y + c = 0
%   y in K*
%
% where K is a cone composed of zero variables (z), nonnegative orthant
% (l), box cone (bl, bu), second-order (Lorentz) cone (q), and cone of
% positive semidefinite matrices (s); and K* is the dual cone of K.
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
% Lorentz cone, rotated Lorentz cone, and PSD cone can be shifted by a
% lower bound (cb).

opts = obj.opts;

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
Nx.l = (obj.getdimc(Kx,'l'));
Nx.q = (obj.getdimc(Kx,'q'));
Nx.s = (obj.getdimc(Kx,'s'));

% number of constraints per cone type
Na.l = (obj.getdimc(Kc,'l'));
Na.q = (obj.getdimc(Kc,'q'));
Na.s = (obj.getdimc(Kc,'s'));

% rewrite linear and state constraints, that is,
% 
%   lba <= al*x <= uba
%   lbx <=  xl  <= ubx
% 
% into
%
%   -al*x + sba =  0
%   -ac*x + sca = -cba
%     -xl + sbx =  0
%     -xc + scx = -cbx
%
% where sba, sbx are box constrained slack variables 
% and sca, scx are elements of Kc and Kx, respectively.
% To that extent we introduce a new slack variable t with constraint
%
%     0 + t = 1
%
% such that (t,sba,sbx) in Kb = {(t,s) | t*[lba;lbx] <= s <= t*[uba;ubx]}.
Nx_c = n - Nx.l;
Na_c = m - Na.l;

A = [sparse(1,n); -a; -speye(n)];
b = [1; sparse(Na.l,1); -cba; sparse(Nx.l,1); -cbx];
c = g;
P = h;

% set cone for SCS
K.bl = [lba; lbx];
K.bu = [uba; ubx];
K.q = [Na.q Nx.q];
K.s = [Na.s Nx.s];

obj.cone = casadi.Function('K',struct2cell(obj.args_in),struct2cell(K),fieldnames(obj.args_in),fieldnames(K));

% reorder slack variables to (b,q,s)
Ia = [ones(1,Na.l) 2*ones(1,sum(Na.q)) 3*ones(1,sum(Na.s.^2))];
Ix = [ones(1,Nx.l) 2*ones(1,sum(Nx.q)) 3*ones(1,sum(Nx.s.^2))];

[~,idx] = sort([0 Ia Ix]);

% SCS-style vectorization of PSD cone
Ac = mat2cell(A(idx,:),[1+Na.l+Nx.l sum(K.q) sum(K.s.^2)],n);
bc = mat2cell(b(idx)  ,[1+Na.l+Nx.l sum(K.q) sum(K.s.^2)],1);
% vector-based vectorization
Ac_s = cellfun(@obj.sdp_vec, mat2cell(Ac{3},K.s.^2,n), 'UniformOutput', false);
bc_s = cellfun(@obj.sdp_vec, mat2cell(bc{3},K.s.^2,1), 'UniformOutput', false);
% build SCS structures
Ascs = vertcat(Ac{1},Ac{2},Ac_s{:});
bscs = vertcat(bc{1},bc{2},bc_s{:});

% return SCS structures P,A,b,c
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),{P Ascs bscs c},fieldnames(obj.args_in),{'P' 'A' 'b' 'c'});

% parse SCS solution (X,Y,S) into (x,cost,lam_a,lam_x)
X = casadi.MX.sym('x',size(c));
Y = casadi.MX.sym('y',size(bscs));
S = casadi.MX.sym('s',size(bscs));

% dual variables corresponding to constraints and variables
Yc = mat2cell(Y,[1 Na.l Nx.l sum(Na.q) sum(Nx.q) sum(Na.s.*(Na.s+1)/2) sum(Nx.s.*(Nx.s+1)/2)],1);
% multipliers for box constraints
lam_a_l = -Yc{2};
lam_x_l = -Yc{3};
% multipliers for Lorentz cone constraints
lam_a_q = -Yc{4};
lam_x_q = -Yc{5};
% de-vectorize PSD dual variables
Yc_a_s = cellfun(@obj.sdp_mat, mat2cell(Yc{6},Na.s.*(Na.s+1)/2,1), 'UniformOutput', false);
Yc_x_s = cellfun(@obj.sdp_mat, mat2cell(Yc{7},Nx.s.*(Nx.s+1)/2,1), 'UniformOutput', false);
% multipliers for PSD constraints
lam_a_s = -vertcat(Yc_a_s{:});
lam_x_s = -vertcat(Yc_x_s{:});
% build conic multipliers
lam_a = [lam_a_l; lam_a_q; lam_a_s];
lam_x = [lam_x_l; lam_x_q; lam_x_s];

obj.ghan = casadi.Function('g',[struct2cell(obj.args_in)' {X Y S}],{X (X'*(P/2)*X + c'*X) lam_a lam_x},[fieldnames(obj.args_in)' {'X' 'Y' 'S'}],obj.names_out);

end
