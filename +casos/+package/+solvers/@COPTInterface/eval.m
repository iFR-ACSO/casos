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
if ~isfield(opts,'Logging'), opts.Logging = 0; end       % Disable output by default

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
K.f = length(If);                   % Free variables
K.l = K.l - nnz(J) - length(If);    % Positive orthants 
                                    % (Adjust for removed variables)

% determine if the problem is an LP
isLP = all( cellfun(@(x) isempty(x) || all(x == 0), {K.q, K.s, K.r}) );
    
if isLP
    % LP formulation
    problem.objsen = 'Minimize';
    problem.A      = A;
    problem.obj    = c;
    problem.lb     = [-inf(K.f,1); zeros(K.l,1)];
    problem.ub     = [+inf(K.f,1); +inf(K.l,1)];
    problem.sense  = repmat('E', size(b));
    problem.rhs    = b;
else
    % general conic formulation
    problem.conedata = struct(   ...
            'objsen',   'min',   ...
            'objcon',   0,       ...
            'K',        K,       ...
            'c',        c,       ...
            'A',        A,       ...
            'b',        b        ...
    );
end

% solve the optimization problem
copt_sol = copt_solve(problem, opts);

% save the info
obj.info = struct( ...
        'status',       copt_sol.status,      ...
        'simplexiter',  copt_sol.simplexiter, ...
        'barrieriter',  copt_sol.barrieriter, ...
        'solvingtime',  copt_sol.solvingtime  ...
);

% error handling for infeasibility
if strcmp(obj.info.status, 'infeasible')
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Conic problem is primal infeasible.')
    x_copt = zeros(size(A,2),1);
    s_copt = zeros(size(A,2),1);
elseif strcmp(obj.info.status, 'unbounded') 
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Conic problem is primal unbounded.')
    x_copt = zeros(size(A,2),1);
    s_copt = zeros(size(A,2),1);
elseif strcmp(obj.info.status, 'numerical') 
   obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
   assert(~obj.opts.error_on_fail,'Conic problem has encountered numerical problems.')
   x_copt = zeros(size(A,2),1);
   s_copt = zeros(size(A,2),1);
else
    % handle success case
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;

    % safely get PSD parts if they exist, otherwise empty
    psdx   = [];
    psdrc  = [];    
    if isfield(copt_sol, 'psdx'),  psdx  = copt_sol.psdx;  end
    if isfield(copt_sol, 'psdrc'), psdrc = copt_sol.psdrc; end

    % assign full solution
    x_copt = [copt_sol.x; psdx];
    s_copt = [copt_sol.rc; psdrc];
   
end

% assign full solution
x = sparse(idx,1,x_copt,length(J),1);
s = sparse(idx,1,s_copt,length(J),1);

% parse solution
argout = call(obj.ghan,[argin {x s}]);

end

