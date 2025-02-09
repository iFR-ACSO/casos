function argout = eval(obj,argin)
% Call COPT interface.

% pre-process bounds
lba = sparse(argin{4});
uba = sparse(argin{5});
cba = sparse(argin{6});
lbx = sparse(argin{7});
ubx = sparse(argin{8});
cbx = sparse(argin{9});

% dimensions of original problem
nl = length(lbx);
nc = length(cbx);
ml = length(lba);
mc = length(cba);

% detect infinite lower variable bounds
If = find(isinf(lbx));
% prepare for infinite lower variable bounds
argin{7}(If) = 0;

% evaluate problem structure
prob = call(obj.fhan,argin);

% to double
A = sparse(prob{1});
b = sparse(prob{2});
c = sparse(prob{3});
% cone
K = obj.cone;

% options to COPT
opts = obj.opts.copt;
% disable output by default
if ~isfield(opts,'fid'), opts.fid = 0; end

% reorder decision variables
idx = [If' setdiff(1:length(c),If)];

% remove trivial constraints
I = false(size(b));
J = false(size(c));

% detect equality constraints
Ila = find(lba == uba);
% remove lower bound constraints
I(Ila) = true;
% remove slack variables (sua,sla)
J(nl+[Ila; ml+Ila]) = true;

% detect constant variables
Ilx = find(lbx == ubx);
% remove slack variable sux
% (this leaves zl(i) = 0)
J(nl+2*ml+Ilx) = true;

% detect infinite bounds
Iba = find(isinf([uba;lba]));
Ibx = find(isinf(ubx));
% remove infinite bound constraints
I(Iba) = true;
I(2*ml+mc+Ibx) = true;
% remove slack variables (sua,sla,sux)
J(nl+[Iba; 2*ml+Ibx]) = true;

% purge constraints
A(I,:) = [];
b(I)   = [];

% purge variables
idx(J) = [];

A = A(:,idx);
c = c(idx);

% modify cone
K.f = length(If);
K.l = K.l - nnz(J) - length(If);

% try COPT ----------------------------

[Alin,Abar] = separate(A, size(A,1), [K.f sum(K.s.^2)]);
Abar_vec = sdp_vec2(Abar, K.s, 1, 2);

[clin,cbar] = separate(c, [K.f sum(K.s.^2)], 1);
cbar_vec = sdp_vec2(cbar, K.s, 1);

problem.conedata.objsen = 'min';
problem.conedata.objcon = 0;
problem.conedata.K.f    = K.f;  % Free variables
problem.conedata.K.l    = K.l;  % Positive orthants
problem.conedata.K.q    = K.q;  % Quadratic cones
problem.conedata.K.r    = K.r;  % Rotated cones
problem.conedata.K.s    = K.s;  % SDP cones
problem.conedata.c      = full([clin; cbar_vec]);
problem.conedata.A      = [Alin, Abar_vec];
problem.conedata.b      = full(b);

% Set parameter
parameter.TimeLimit = 60;
parameter.Logging = 0;
solution = copt_solve(problem, parameter);

% save the info
obj.info.status         = solution.status;
obj.info.simplexiter    = solution.simplexiter;
obj.info.barrieriter    = solution.barrieriter;
obj.info.solvingtime    = solution.solvingtime;
obj.info.rowmap         = solution.rowmap;

if strcmp(obj.info.status, 'infeasible')
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Conic problem is primal infeasible.')
    x_copt = zeros(length(c),1);
    y_copt = zeros(size(A,1),1);
else
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
    x_copt = [solution.x];
    len_size = K.s.*(K.s+1)/2;
    for i=1:length(K.s)
        temp = unflatten_sdp(...
            solution.psdx((1+sum(len_size(1:(i-1)))):sum(len_size(1:i))), K.s(i));
        x_copt = [x_copt; temp(:)];
    end
    y_copt = solution.psdpi;
end

% -------------------------------------
% call SeDuMi
%tic
%[x_,y_,obj.info] = sedumi(A,b,full(c),K,opts);
%toc

% assign full solution
x = sparse(idx,1,x_copt,length(J),1);
y = sparse(find(~I),1,y_copt,length(I),1);

% parse solution
argout = call(obj.ghan,[argin {x y}]);

end

function varargout = separate(A,varargin)
% Separate array into subarrays.

    varargout = mat2cell(A,varargin{:});
end

function [v,i,j,k,l] = sdp_vec2(M,Ks,scale,dim)
% Index-based lower-triangular vectorization for semi-definite matrices.
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
% corresponds to stacking the lower-triangular elements of Mij column-wise.
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

assert(size(M,dim) == sum(Nq), 'Number and size of block entries in M must correspond to dimensions Ks.')

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
    V = zeros(sz);
    
    % Assign values to V using linear indexing
    ind = sub2ind(sz, subIv{:});  
    V(ind) = val;  % Fill with given values
    
    % Return block matrix V
    v = V;
end

end

function X = unflatten_sdp(x_flat, n)
    % Convert a flattened symmetric matrix representation back to a full matrix.
    %
    % INPUTS:
    %   x_flat - Vector of independent elements of a symmetric n x n matrix
    %   n      - Dimension of the original matrix
    %
    % OUTPUT:
    %   X      - Reconstructed n x n symmetric matrix

    % Initialize empty symmetric matrix
    X = zeros(n, n);

    % Index counter for flattened vector
    idx = 1;
    
    for i = 1:n
        for j = i:n
            X(i, j) = x_flat(idx);
            X(j, i) = x_flat(idx); % Use symmetry
            idx = idx + 1;
        end
    end
end