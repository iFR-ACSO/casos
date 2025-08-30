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

% reorder due to PSD(1)
[Mx, new_linx] = reordering_map(n, Nx);     % due to reordering on Kx
[Mc, new_linc] = reordering_map(m, Na);     % due to reordering on Kc

a = Mc*a*Mx'; % perform permuation
g = Mx*g;     % perform permutaion

% extend lbx and ubx if PSD(1) existed in Kx
lbx = [lbx; zeros(new_linx)];
ubx = [ubx;   inf(new_linx)];

% extend lba and uba if PSD(1) existed in Kc
lba = [lba; zeros(new_linc)];
uba = [uba;   inf(new_linc)];

% update cbx due to variable reordering
[cbxq, cbxr, cbxs]= separate(cbx, [sum(Nx.q); sum(Nx.r); sum(Nx.s.^2)]);
cbx_split = cell(length(Nx.s), 1);
[cbx_split{:}] = separate(cbxs, Nx.s.^2); 
cbx = [cbxq; cbxr; vertcat(cbx_split{Nx.s~=1})];

% handle PSD cones of size 1 (change cone sizes)
Nx.l = Nx.l + sum(Nx.s == 1);
Nx.s = Nx.s(Nx.s > 1);

Na.l = Na.l + sum(Na.s == 1);
Na.s = Na.s(Na.s > 1);

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

% set cone for COPT
K.l = 2*(Na.l + Nx.l);
K.q = [Na.q Nx.q];
K.r = [Na.r Nx.r];
K.s = [Na.s Nx.s];

% Check if any entry in K.s is equal to 1, and throw an error if so
assert(~any(K.s == 1), 'COPT does not support PSD cones of size 1. Use a linear cone instead.');

obj.cone = K;

% reorder decision variables (z,s) to (l,q,r,s)
% current order: (z,s) = [zl zc sua sla sca sux]
Iz = [ones(1,Nx.l) 2*ones(1,sum(Nx.q)) 3*ones(1,sum(Nx.r)) 4*ones(1,sum(Nx.s.^2))];
Ia = [ones(1,2*Na.l) 2*ones(1,sum(Na.q)) 3*ones(1,sum(Na.r)) 4*ones(1,sum(Na.s.^2))];
Ix =  ones(1,Nx.l);

% new order: [zl sua sla sux | zc_q sca_q zc_r sca_r zc_s sca_s] 
% (see https://de.mathworks.com/help/matlab/ref/sort.html#bt8nojg-1-I)
[~,idx] = sort([Iz Ia Ix]);

% reorder matrices
A = A(:,idx);
c = c(idx);

% set specific formulation for COPT
% separate A and c based on K
[Alin,Abar] = separate(A, size(A,1), [(K.l), sum(K.s.^2)]);

% vectorize the upper-triangular elements of the SDP blocks
Abar_vec = obj.sdp_vec(Abar, K.s, 1, 2, true);

% separate c into linear part (clin) and semidefinite part (cbar)
[clin,cbar] = separate(c, [(K.l), sum(K.s.^2)], 1);

% vectorize the upper-triangular elements of the SDP blocks
cbar_vec = obj.sdp_vec(cbar, K.s, 1, 1, true);

% Concatenate full COPT vectors
c_copt      = full([clin; cbar_vec]);   % vector for objective
A_copt      = [Alin, Abar_vec];         % constraint matrix
b_copt      = full(b);                  % constraint right-hand side

% return COPT structures A,b,c
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),{A_copt b_copt c_copt},fieldnames(obj.args_in),{'A' 'b' 'c'});

% define CasADi symbolic variables for the primal problem
% decision variables
X = casadi.MX.sym('x',size(c_copt));
% slack/dual variables
S = casadi.MX.sym('s',size(c_copt));



% separate X into linear and PSD components
[Xlin, Xpsd] = separate(X, [size(Alin,2) size(Abar_vec,2)], 1);
% reshape the PSD component into proper SDP matrix form
Xpsd = obj.sdp_mat(Xpsd,  [Nx.s Na.s], 1, 1, false);    
Xfull = [Xlin; Xpsd];   % combine linear and PSD parts

% separate S into linear and PSD components
[Slin, Spsd] = separate(S, [size(Alin,2) size(Abar_vec,2)], 1);
% reshape the PSD component into proper SDP matrix form
Spsd = obj.sdp_mat(Spsd,  [Nx.s Na.s], 1, 1, false);  
Sfull = [Slin; Spsd];   % combine linear and PSD parts

% primal decision variables in the NEW order
Xc = mat2cell(Xfull,[Nx.l Na.l Na.l Nx.l sum(Nx.q) sum(Na.q) sum(Nx.r) sum(Na.r) sum(Nx.s.^2) sum(Na.s.^2)],1);
x = [Xc{1}; Xc{5}; Xc{7}; Xc{9}] + zb;

% dual variables corresponding to constraints and variables
% S in the ORIGINAL order of (z,s), 
% that is, S = [y_lbx y_cbx y_uba y_lba y_cba y_ubx]
Sc = mat2cell(Sfull,[Nx.l Nx_c Na.l Na.l Na_c Nx.l],1);

% multipliers for interval constraints
lam_a_l = Sc{3} - Sc{4};    % y_uba - y_lba
lam_x_l = Sc{6} - Sc{1};    % y_ubx - y_lbx

% multipliers for cone constraints
lam_a_c = -Sc{5};  %  -y_cba
lam_x_c = -Sc{2};  %  -y_cbx

obj.ghan = casadi.Function('g',[struct2cell(obj.args_in)' {X S}],{Mx'*x (g'*x) Mc'*[lam_a_l; lam_a_c] Mx'*[lam_x_l; lam_x_c]},[fieldnames(obj.args_in)' {'X' 'S'}],obj.names_out);

end


function varargout = separate(A,varargin)
% Separate array into subarrays.
    varargout = mat2cell(A,varargin{:});
end

function [M, new_lin] = reordering_map(size, N)
% Moves scalar PSD cones to the linear part and keeps other cones in order.
%
% Inputs:
%   size - total number of variables
%   N    - structure with fields l, q, r, s (number of variables per cone)
%
% Outputs:
%   M        - sparse permutation matrix
%   new_lin  - number of scalar PSD cones moved to linear part

    % start with identity mapping
    M = speye(size);

    % sizes of each group
    sz_l = sum(N.l);
    sz_q = sum(N.q);
    sz_r = sum(N.r);
    sz_s = sum(N.s.^2);

    % separate the mapping into each component
    [Ml, Mq, Mr, Ms]= separate(M, [sz_l; sz_q; sz_r; sz_s]);
    
    % split each PSD cone
    Ms_split = cell(length(N.s), 1);
    [Ms_split{:}] = separate(Ms, N.s.^2); 
    
    % indices of scalar PSD cones
    psd1_idx = (N.s == 1);
    
    % build the mapping
    M = vertcat(Ml,                             ...
                vertcat(Ms_split{psd1_idx}),    ...    % all size-1 PSD cones
                Mq,                             ...
                Mr,                             ...
                vertcat(Ms_split{~psd1_idx}));         % all larger PSD cones
    
    % number of scalar PSD cones moved to linear
    new_lin = nnz(psd1_idx);
end
