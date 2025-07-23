function argout = eval(obj,argin)
% Evaluate the optimization problem using the COPT interface.

% Pre-process bounds (convert to sparse format)
lba = sparse(argin{4});
uba = sparse(argin{5});
cba = sparse(argin{6});
lbx = sparse(argin{7});
ubx = sparse(argin{8});
cbx = sparse(argin{9});

% Determine original problem dimensions
nl = length(lbx);
nc = length(cbx);
ml = length(lba);
mc = length(cba);

% Detect and handle infinite lower variable bounds
If = find(isinf(lbx));
argin{7}(If) = 0;       % Replace infinities with zero for processing

% Evaluate the problem structure
prob = call(obj.fhan,argin);

% Convert problem data to sparse double format
A = sparse(prob{1});
b = sparse(prob{2});
c = sparse(prob{3});

% Extract cone structure
K = obj.cone;

% Extract options for COPT solver
opts = obj.opts.copt;

% Set default solver parameters if not already specified
if ~isfield(opts,'TimeLimit'), opts.TimeLimit = 60; end  % Default time limit
if ~isfield(opts,'Logging'), opts.Logging = 0; end       % Disable output by default

% Save size of original decision variable vector (vectorized)
len_orig = sum([K.l])+sum(K.s.*K.s);

% Reorder decision variables
idx = [If' setdiff(1:length(c),If)];

% Initialize constraint and variable masks
I = false(size(b));     % Constraints to remove
J = false(size(c));     % Variables to remove

% Detect and remove equality constraints
Ila = find(lba == uba); % Find equalities
I(Ila) = true;          % remove lower bound constraints

% Remove slack variables (sua, sla)
J(nl+[Ila; ml+Ila]) = true;

% Detect and remove constant variables
Ilx = find(lbx == ubx); % Find constant variables
J(nl+2*ml+Ilx) = true;  % remove slack variable sux

% Detect and remove infinite bounds
Iba = find(isinf([uba;lba]));
Ibx = find(isinf(ubx));

% Remove infinite bound constraints
I(Iba) = true;
I(2*ml+mc+Ibx) = true;

% Remove slack variables (sua,sla,sux)
J(nl+[Iba; 2*ml+Ibx]) = true;
 
% Purge constraints
A(I,:) = [];
b(I)   = [];

% Purge variables
idx(J) = [];
A = A(:,idx);
c = c(idx);

% Modify cone structure
K.f = length(If);
K.l = K.l - nnz(J) - length(If);    % Adjust for removed variables

% Setup problem structure
problem.conedata = struct( ...
        'objsen', 'min', ...
        'objcon', 0, ...
        'K', struct( ...
                    'f', K.f, ...   % Free variables
                    'l', K.l, ...   % Positive orthants
                    'q', K.q, ...   % Quadratic cones
                    'r', K.r, ...   % Rotated cones
                    's', K.s  ...   % SDP cones
                    ), ...
        'c', full(c), ...
        'A', A, ...
        'b', full(b) ...
);

% Solve the optimization problem
solution = copt_solve(problem, opts);

% save the info
obj.info = struct( ...
        'status', solution.status, ...
        'simplexiter', solution.simplexiter, ...
        'barrieriter', solution.barrieriter, ...
        'solvingtime', solution.solvingtime, ...
        'rowmap', solution.rowmap ...
);


idx = [If' setdiff(1:len_orig,If)];  % reorder decision variables
J = false([len_orig,1]);             % remove trivial constraints
Ila = find(lba == uba);              % detect equality constraints
J(nl+[Ila; ml+Ila]) = true;          % remove slack variables (sua,sla)
Ilx = find(lbx == ubx);              % detect constant variables
J(nl+2*ml+Ilx) = true;               % remove slack variable sux
Iba = find(isinf([uba;lba]));        % detect infinite bounds
Ibx = find(isinf(ubx));
I(Iba) = true;                       % remove infinite bound constraints
I(2*ml+mc+Ibx) = true;
J(nl+[Iba; 2*ml+Ibx]) = true;        % remove slack variables (sua,sla,sux)
idx(J) = [];                         % purge variables

% Error handling for infeasibility
if strcmp(obj.info.status, 'infeasible')
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Conic problem is primal infeasible.')
    x_copt = zeros(length(idx),1);
    y_copt = zeros(size(A,1),1);
elseif strcmp(obj.info.status, 'unbounded') 
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Conic problem is primal unbounded.')
    x_copt = zeros(length(idx),1);
    y_copt = zeros(size(A,1),1);
elseif strcmp(obj.info.status, 'numerical') 
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
    assert(~obj.opts.error_on_fail,'Conic problem has encountered numerical problems.')
    x_copt = zeros(length(idx),1);
    y_copt = zeros(size(A,1),1);
else
    % Success case
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
    
    % Get from solution the primal and dual variables
    x_temp = [solution.x];
    y_temp = [solution.pi];

    % Check if fields exist and are non-empty
    if isfield(solution, 'psdx') && ~isempty(solution.psdx)
        x_temp = [x_temp; solution.psdx];
    end
    
    if isfield(solution, 'psdpi') && ~isempty(solution.psdpi)
        y_temp = solution.psdpi;
    end
	
    % Extract primal variables
    num_x = sum([K.f, K.l]); 
    x_copt = zeros(num_x + sum(K.s.^2), 1);  % Preallocate to final size
        
    x_copt(1:num_x) = x_temp(1:num_x);
    psdx = x_temp(num_x + 1:end);


    % Unflatten SDP solution into vector form
    len_size = K.s .* (K.s + 1) / 2;        % Size of each SDP block in vector form
    offset = num_x;                         % Track position in x_copt

    for i = 1:length(K.s)
        temp = unflatten_sdp(psdx(sum(len_size(1:i-1)) + 1 : sum(len_size(1:i))), K.s(i));
        x_copt(offset + (1:numel(temp))) = temp(:);
        offset = offset + numel(temp);
    end
    y_copt = y_temp;

end

% assign full solution
x = sparse(idx,1,x_copt,length(J),1);
y = sparse(find(~I),1,y_copt,length(I),1);

% parse solution
argout = call(obj.ghan,[argin {x y}]);

end


function X = unflatten_sdp(x_flat, n)
    % Convert a flattened symmetric matrix representation back to a full matrix.
    % INPUTS:
    %   x_flat - Vector of independent elements of a symmetric n x n matrix
    %   n      - Dimension of the original matrix
    % OUTPUT:
    %   X      - Reconstructed n x n symmetric matrix

    % TODO: Fix bug with more efficient code below
    % [row, col] = meshgrid(1:n, 1:n);
    % mask = row <= col;
    % idx = 1;
    % X(mask) = x_flat(idx:idx+sum(mask(:))-1);
    % X = X + X' - diag(diag(X));

    % OLD:
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
