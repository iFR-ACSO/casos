function buildproblem(obj)
% Construct COPT problem description (A,b,c,K) from problem structure.
%
% COPT solves semidefinite problems in the form
%
%   min c'*x st. A*x = b and x in K
%
% with dual variables satisfying 
%
%   c - A'*y in K*
%
% where K is a cone composed of free variables (f), nonnegative orthant
% (l), Lorentz cone (q), rotated Lorentz cone (r), and cone of positive
% semidefinite matrices (s); and K* is the dual cone of K.
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
Nx.l = (obj.getdimc(Kx,'lin'));
Nx.q = (obj.getdimc(Kx,'lor'));
Nx.r = (obj.getdimc(Kx,'rot'));
Nx.s = (obj.getdimc(Kx,'psd'));

% number of constraints per cone type
Na.l = (obj.getdimc(Kc,'lin'));
Na.q = (obj.getdimc(Kc,'lor'));
Na.r = (obj.getdimc(Kc,'rot'));
Na.s = (obj.getdimc(Kc,'psd'));

% define xl = zl + lbx, xc = zc + cbx, 
% where zl is a nonnegative variable and zc is element of Kx
zb = [lbx; cbx];
% rewrite linear and state constraints, that is,
%
%   lba <= al*x <= uba, ac*x in Kc
%   lbx <=  xl  <= ubx,  xc  in Kx
%
% with a = [al; ac], x = (xl, xc)
% into equality constraints
%
%   al*x + sua = uba <-> al*z + sua = uba - al*zb
%   al*x - sla = lba <-> al*z - sla = lba - al*zb
%   ac*x - sca = cba <-> ac*z - sca = cba - ac*zb
%     xl + sux = ubx <->   zl + sux = ubx -   lbx
%
% where sla, sua, and sux are nonnegative slack variables 
% and sca is element of Kc
Nx_c = n - Nx.l;
Na_c = m - Na.l;

% [al; al; ac]
A0 = [a(1:Na.l,:); a];

% min c'*(z,s) s.t. A*(z,s) = b
A = [ [A0; speye(Nx.l,n)] blkdiag(+speye(Na.l), -speye(m), +speye(Nx.l)) ];
b = [ [uba; lba; cba] - A0*zb; ubx - lbx ];
c = [g; sparse(m+Na.l+Nx.l,1)];

% set cone for SeDuMi
K.l = 2*(Na.l + Nx.l);
K.q = [Na.q Nx.q];
K.r = [Na.r Nx.r];
K.s = [Na.s Nx.s];

obj.cone = K;

% reorder decision variables (z,s) to (l,q,r,s)
% current order: (z,s) = [zl zc sua sla sca sux]
Iz = [ones(1,Nx.l) 2*ones(1,sum(Nx.q)) 3*ones(1,sum(Nx.r)) 4*ones(1,sum(Nx.s.^2))];
Ia = [ones(1,2*Na.l) 2*ones(1,sum(Na.q)) 3*ones(1,sum(Na.r)) 4*ones(1,sum(Na.s.^2))];
Ix =  ones(1,Nx.l);

% new order: [zl sua sla sux | zc_q sca_q zc_r sca_r zc_s sca_s] 
% (see https://de.mathworks.com/help/matlab/ref/sort.html#bt8nojg-1-I)
[~,idx] = sort([Iz Ia Ix]);

% return COPT structures A,b,c
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),{A(:,idx) b c(idx)},fieldnames(obj.args_in),{'A' 'b' 'c'});

% parse COPT solution (X,Y) into (x,cost,lam_a,lam_x)
X = casadi.MX.sym('x',size(c));
Y = casadi.MX.sym('y',size(b));
S = c - A'*Y;

% primal decision variables in the NEW order
Xc = mat2cell(X,[Nx.l Na.l Na.l Nx.l sum(Nx.q) sum(Na.q) sum(Nx.r) sum(Na.r) sum(Nx.s.^2) sum(Na.s.^2)],1);
x = [Xc{1}; Xc{5}; Xc{7}; Xc{9}] + zb;

% dual variables corresponding to constraints and variables
% S in the ORIGINAL order of (z,s), 
% that is, S = [y_lbx y_cbx y_uba y_lba y_cba y_ubx]
Sc = mat2cell(S,[Nx.l Nx_c Na.l Na.l Na_c Nx.l],1);
% multipliers for interval constraints
lam_a_l = Sc{3} - Sc{4};    % y_uba - y_lba
lam_x_l = Sc{6} - Sc{1};    % y_ubx - y_lbx
% multipliers for cone constraints
lam_a_c = -Sc{5};  %  -y_cba
lam_x_c = -Sc{2};  %  -y_cbx

obj.ghan = casadi.Function('g',[struct2cell(obj.args_in)' {X Y}],{x (g'*x) [lam_a_l; lam_a_c] [lam_x_l; lam_x_c]},[fieldnames(obj.args_in)' {'X' 'Y'}],obj.names_out);

end
