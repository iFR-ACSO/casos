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

% Check if any entry in K.s is equal to 1, and throw an error if so
assert(~any(K.s == 1), 'COPT does not support PSD cones of size 1. Use a linear cone instead.');

obj.cone = K;

% reorder decision variables (z,s) to (l,q,r,s)
% current order: (z,s) = [zl zc sua sla sca sux]
Iz = [ones(1,Nx.l) 2*ones(1,sum(Nx.q)) 3*ones(1,sum(Nx.r)) 4*ones(1,sum(Nx.s.^2))];
Ia = [ones(1,2*Na.l) 2*ones(1,sum(Na.q)) 3*ones(1,sum(Na.r)) 4*ones(1,sum(Na.s.^2))];
%Iz = [ones(1,Nx.l) 2*ones(1,sum(Nx.q)) 3*ones(1,sum(Nx.r)) 4*ones(1,sum(Nx.s.*(Nx.s+1)./2))];
%Ia = [ones(1,2*Na.l) 2*ones(1,sum(Na.q)) 3*ones(1,sum(Na.r)) 4*ones(1,sum(Na.s.*(Na.s+1)./2))];
Ix =  ones(1,Nx.l);

% new order: [zl sua sla sux | zc_q sca_q zc_r sca_r zc_s sca_s] 
% (see https://de.mathworks.com/help/matlab/ref/sort.html#bt8nojg-1-I)
[~,idx] = sort([Iz Ia Ix]);

% reorder matrices
A = A(:,idx);
c = c(idx);

% set specific formulation for COPT
% Separate A and c based on K
[Alin,Abar] = separate(A, size(A,1), [(K.l), sum(K.s.^2)]);
Abar_vec = sdp_vec_upper(Abar, K.s, 1, 2);

[clin,cbar] = separate(c, [(K.l), sum(K.s.^2)], 1);
cbar_vec = sdp_vec_upper(cbar, K.s, 1);
c_copt      = full([clin; cbar_vec]);
A_copt      = [Alin, Abar_vec];
b_copt      = full(b);
%c_copt      = c;
%A_copt      = A;
%b_copt      = b;

% return COPT structures A,b,c
%obj.fhan = casadi.Function('f',struct2cell(obj.args_in),{A_copt(:,idx) b_copt c_copt(idx)},fieldnames(obj.args_in),{'A' 'b' 'c'});
obj.fhan = casadi.Function('f',struct2cell(obj.args_in),{A_copt b_copt c_copt},fieldnames(obj.args_in),{'A' 'b' 'c'});

% parse COPT solution (X,Y) into (x,cost,lam_a,lam_x)
%X_copt = casadi.MX.sym('x',size(c_copt));
%Y_copt = casadi.MX.sym('y',size(b_copt));

X_copt = casadi.MX.sym('x',size(c));
Y_copt = casadi.MX.sym('y',size(b));

% Unflatten SDP solution into vector form
%len_size = K.s.*(K.s+1)/2;  % Size of each SDP block in vector form
%X = X_copt(1:K.l);
%for i=1:length(K.s)
%    temp = unflatten_sdp(...
%           X_copt((1+sum(len_size(1:(i-1)))):sum(len_size(1:i))), K.s(i));
%    X = [X; temp(:)];
%end
X = X_copt;
Y = Y_copt;
S = c - A'*Y;

%[zblin,zbbar] = separate(zb, [(Nx.l), sum(K.s.^2)], 1);
%zbbar_vec = sdp_vec(zbbar, K.s, 1);
%zb      = full([zblin; zbbar_vec]);

% primal decision variables in the NEW order
Xc = mat2cell(X,[Nx.l Na.l Na.l Nx.l sum(Nx.q) sum(Na.q) sum(Nx.r) sum(Na.r) sum(Nx.s.^2) sum(Na.s.^2)],1);
%Xc = mat2cell(X_copt,[Nx.l Na.l Na.l Nx.l sum(Nx.q) sum(Na.q) sum(Nx.r) sum(Na.r) sum(Nx.s.*(Nx.s+1)./2) sum(Na.s.*(Na.s+1)./2)],1);
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

obj.ghan = casadi.Function('g',[struct2cell(obj.args_in)' {X_copt Y_copt}],{x (g'*x) [lam_a_l; lam_a_c] [lam_x_l; lam_x_c]},[fieldnames(obj.args_in)' {'X' 'Y'}],obj.names_out);

end


function varargout = separate(A,varargin)
% Separate array into subarrays.

    varargout = mat2cell(A,varargin{:});
end

function [v,i,j,k,l] = sdp_vec_upper(M,Ks,scale,dim)
    % Index-based upper-triangular vectorization for semi-definite matrices.
    %
    % This function takes a matrix
    %
    %       | m11 ... m1q |
    %   M = |  :       :  |
    %       | mp1 ... mpq |
    %
    % where each block entry mij is either a column or row vector satisfying
    % mij = Mij(:) for some N-by-N matrix Mij, and computes the matrix
    %
    %       | v11 ... v1q |
    %   V = |  :       :  |
    %       | vp1 ... vpq |
    %
    % where each block entry vij is a N*(N+1)/2-by-1 (column/row) vector that
    % corresponds to stacking the upper-triangular elements of Mij column-wise.
    %
    % Syntax #1:
    %
    %   V = sdp_vec(M,Ks,scale,dim)
    %
    % Returns the block matrix V as described above with the following
    % parameters:
    %
    % - dim:    Whether the blocks of M and V are treated as columns (dim = 1) 
    %           or row vectors (dim = 2); default: dim = 1 if M is a matrix.
    % - scale:  Number to scale the off-diagonal terms of each Mij by; 
    %           default: scale = sqrt(2).
    % - Ks:     Dimensions Nij of the matrices Mij; if dim = 1, then Ks is a
    %           p-by-1 vector satisfying Nij = K(i); otherwise, Ks is a q-by-1
    %           vector satisfying Nij = K(j) for all (i,j) in {1...p}x{1...q}.
    %
    % Syntax #2:
    %
    %   [val,i,j,k,l] = sdp_vec(M,...)
    %
    % Returns the nonzero elements of V along with vectors of indices (i,j) and
    % (k,l) corresponding to the lower-triangular elements Mij(k,l). Parameters
    % apply as above.
    
    if nargin > 3
        % dimension provided, nothing to do
    elseif isrow(M)
        % row vector blocks only
        dim = 2;
    else
        % column-only or default for matrix
        dim = 1;
    end
    
    if nargin < 3 || isempty(scale)
        % default scaling for SCS
        scale = sqrt(2);
    end
    
    % ensure matrix dimensions are a row vector
    s = reshape(Ks,1,[]);
    % number of elements in each matrix block
    Nq = s.^2;
    
    assert(size(M,dim) == sum(Nq), ...
        'Number and size of block entries in M must correspond to dimensions Ks.')
    
    % get nonzero subindices w.r.t. block matrix M
    subI = cell(1,2);
    [subI{:}] = ind2sub(size(M),1:numel(M));
    
    % select dimension and ensure linear indices are a row vector
    I = reshape(subI{dim},1,[]);
    
    % number of elements before matrix variables Mij
    S = cumsum([0 Nq(1:end-1)]);
    
    % compute off-dimension index of corresponding matrix variable
    J = sum(I > S', 1); % interface is 1-based
    
    % compute linear indices of elements in each matrix
    I0 = I - S(J);
    
    % compute indices (k,l) of elements in matrix Mij, where I0 = N*l + k;
    % as (remainder,quotient) of linear indices divided by matrix size
    l = ceil(I0 ./ s(J)); % col
    k = I0 - s(J).*(l-1); % row
    
    % determine strictly upper triangle
    triu = (k < l);
    
    % remove indices for strictly upper triangle
    J(triu) = [];
    l(triu) = [];
    k(triu) = [];
    
    % linear indices for lower triangular entries
    subItril = {subI{1}(~triu) subI{2}(~triu)};
    Itril = sub2ind(size(M),subItril{:});
    
    % determine strictly lower and upper triangle
    tril = (k > l);
    
    % scale strictly lower triangle
    scaling = ones(size(tril));
    scaling(tril) = scale;
    
    % nonzero elements of lower triangular matrices, scaled
    val = scaling.*reshape(M(Itril),1,length(Itril));
    
    if nargout > 1
        % return subindices (i,j) of matrices Mij
        subij{dim} = J; subij{3-dim} = subItril{3-dim};
        [i,j] = subij{:};
        % return nonzero elements
        v = val;
    else
        % Compute linear indices into each block Vij
        kprime = k - l + 1;  % Number of rows from diagonal
        Iv0 = (l-1).*(s(J)-l/2+1) + kprime;  % Number of lower-triangular elements
        
        % Compute cumulative linear indices
        Nv = s.*(s+1)/2; 
        Sv = cumsum([0 Nv(1:end-1)]);
        Iv = Sv(J) + Iv0;
        
        % Subindices into block matrix V
        subIv{dim} = Iv;  % 1-based indexing in MATLAB
        subIv{3-dim} = subItril{3-dim};  
        
        % Size of block matrix V
        sz = size(M); 
        sz(dim) = sum(Nv);
        
        % Initialize the block matrix V as a full (dense) matrix
        V = casadi.MX.zeros(sz);
        
        % Assign values to V using linear indexing
        ind = sub2ind(sz, subIv{:});  
        V(ind) = val;  % Fill with given values
        
        % Return block matrix V
        v = V;
    end
end
