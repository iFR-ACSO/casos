function buildproblem(obj)
% Construct SeDuMi problem description (A,b,c,K) from problem structure.
%
% SeDuMi solves semidefinite problems in the form
%
%   min c'*x st. A*x = b and x in K
%
% where K is a cone composed of free variables (f), nonnegative orthant
% (l), Lorentz cone (q), rotated Lorentz cone (r), and cone of positive
% semidefinite matrices (s).
%
% The generic conic problem in casos has the form
%
%   min 1/2*x'*h*x + g'*x
%   st. a*x in Ka, x in Kx,
%
% where Ka and Kx are cones composed of
%   - box constraints [lb ub] (l)
%   - Lorentz (quadratic) cone (q)
%   - rotated Lorentz cone (r)
%   - cone of PSD matrices (s)
%
% Lorentz cone, rotated Lorentz cone, and PSD cone can be shifted by a
% lower bound (cb).
% Unlike SeDuMi, no free (unrestricted) components are allowed.

opts = obj.sdpopt;

% symbolic variables
g = obj.args_in.g;
a = obj.args_in.a;
% symbolic bounds
lba = obj.args_in.lba;
uba = obj.args_in.uba;
lbx = obj.args_in.lbx;
ubx = obj.args_in.ubx;
% problem size
[m,n] = size(a);

% obtain cones
Kx = opts.Kx;
Ka = opts.Ka;

% number of variables per cone type
Nx.l = (obj.getdimc(Kx,'l'));
Nx.q = (obj.getdimc(Kx,'q'));
Nx.r = (obj.getdimc(Kx,'r'));
Nx.s = (obj.getdimc(Kx,'s'));
% number of constraints per cone type
Na.l = (obj.getdimc(Ka,'l'));
Na.q = (obj.getdimc(Ka,'q'));
Na.r = (obj.getdimc(Ka,'r'));
Na.s = (obj.getdimc(Ka,'s'));

% rewrite linear and state constraints into
%
%   a*x + sua = uba
%   a*x - sla = lba
%   a*x - sca = cba
%     x + sux = ubx
%     x - slx = lbx
%     x - scx = cbx
%
% where sla, sua, slx, sux are nonnegative slack variables and sca, scx are
% elements of Ka and Kx, respectively.
Nx_c = n - Nx.l;
Na_c = m - Na.l;

A = [
    a(1:Na.l,:)    +eye(Na.l) zeros(Na.l,m) zeros(Na.l,n+Nx.l)
    a           zeros(m,Na.l)       -eye(m) zeros(m,   n+Nx.l)
    eye(Nx.l,n) zeros(Nx.l,m+Na.l)           +eye(Nx.l) zeros(Nx.l,n)
    eye(n)      zeros(n,   m+Na.l)        zeros(n,Nx.l)       -eye(n)
];
b = [uba(1:Na.l); lba; ubx(1:Nx.l); lbx];
c = [g; zeros(m+Na.l+n+Nx.l,1)];

% set cone for SeDuMi
K.f = n;
K.l = 2*(Nx.l + Na.l);
K.q = [Nx.q Na.q];
K.r = [Nx.r Na.r];
K.s = [Nx.s Na.s];

obj.cone = K;

% reorder slack variables to (l,q,r,s)
Ia = [ones(1,2*Na.l) 2*ones(1,sum(Na.q)) 3*ones(1,sum(Na.r)) 4*ones(1,sum(Na.s))];
Ix = [ones(1,2*Nx.l) 2*ones(1,sum(Nx.q)) 3*ones(1,sum(Nx.r)) 4*ones(1,sum(Nx.s))];

[~,idx] = sort([zeros(1,n) Ia Ix]);

% return SeDuMi structures A,b,c
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),{A(:,idx) b c(idx)},fieldnames(obj.args_in),{'A' 'b' 'c'});

% parse SeDuMi solution (X,Y) into (x,cost,lam_a,lam_x)
X = casadi.MX.sym('x',size(c));
Y = casadi.MX.sym('y',size(b));

% dual variables corresponding to constraints and variables
Yc = mat2cell(Y,[Na.l Na.l Na_c Nx.l Nx.l Nx_c],1);
% multipliers for interval constraints
lam_a_l = Yc{1} - Yc{2};
lam_x_l = Yc{4} - Yc{5};
% multipliers for cone constraints
lam_a_c = Yc{3};
lam_x_c = Yc{6};

obj.ghan = casadi.Function('g',[struct2cell(obj.args_in)' {X Y}],{X(1:n) (c'*X) [lam_a_l; lam_a_c] [lam_x_l; lam_x_c]},[fieldnames(obj.args_in)' {'X' 'Y'}],obj.names_out);

end
