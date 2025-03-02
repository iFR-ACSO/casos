function Lz = newton_reduce(Pdegmat,Zdegmat,solver)
% Removes monomials outside half Newton polytope
% Strategy inpired from: Simplification Methods for Sum-of-Squares 
% Programs, Peter Seiler et al.

% list of available solvers for the newton polytope reduction
list_solvers = {'linprog', 'sedumi', 'mosek', 'scs'};
[solver_available, solver_id] = ismember(solver,list_solvers);

if solver_available==0
    error('Specified solver is not available.\n');
end

switch solver_id
    case 1  % options for linprog (Optimization toolbox from Matlab)
        options = optimoptions('linprog');
        options.Display = 'off';   
        lp_solver = @(A,b,c,options) linprog_solver(A, b, c, options);
    case 2  % options for Sedumi
        options.fid = 0;
        lp_solver = @(A,b,c,options) sedumi_lp_solver(A, b, c, options);
    case 3  % options for Mosek
        options = struct();
        lp_solver = @(A,b,c,options) mosek_lp_solver(A, b, c, options);
    case 4 % options for SCS
        options = struct();
        lp_solver = @(A,b,c,options) scs_lp_solver(A, b, c, options);
end

% build LP (part 1)
bfixed = [zeros(size(Pdegmat,1),1); 1];

% trivial removal of monomials
keep_trivial = ismember(Zdegmat*2,Pdegmat,'rows');
keep = true(size(keep_trivial));

% try to go over each possible monomial basis and verify if it belongs to
% the newton polytope by checking for the existance of a hyperplane 
for i = 1:length(keep_trivial)

    if keep_trivial(i) || ~keep(i)
        continue;
    end
    
    % build LP (part 2)
    q = Zdegmat(i,:)*2;
    c = ([-q 1])';
    F_struc = ([bfixed -[Pdegmat -ones(size(Pdegmat,1),1); q -1]]);
    
    A = -F_struc(1:end,2:end);
    b = F_struc(1:end,1);   

    % solve LP
    [x,flag] = lp_solver(A,b,c,options);

    % % in case the LP is not feasible, giving an empty output 'x' 
    if isempty(x)
        x = zeros(length(c),1);
    end

    % If the LP gives an unbounded solution or the only solution is a zero vector
    if (flag > 0 && ([-q 1]*x(:) < 0)) || flag == -3
        a = x(1:end-1);
        b1 = x(end);
        u = 2*Zdegmat*a - b1 > sqrt(eps);
        keep(u) = 0;
    end

end

Lz = keep;
Lz = sparse(Lz);
end


% -------------------------------------------------------------------------
function [x, flag] = mosek_lp_solver(A,b,c,options)
    % MOSEK_LP_SOLVER Solves a linear program using the MOSEK solver
    %   Minimize:   c' * x
    %   Subject to: A * x <= b
    %
    % Inputs:
    %   A        - Constraint matrix (m x n)
    %   b        - Right-hand side vector (m x 1)
    %   c        - Objective function coefficients (n x 1)
    %   options  - Struct with optional solver settings (optional)
    %
    % Outputs:
    %   x        - Solution vector
    %   flag     - Solver status: 
    %               1  -> Solved
    %              -3  -> Infeasible
    %              -2  -> Unbounded
    %               0  -> Other cases (failure)

    % Validate inputs
    if nargin < 3
        error('mosek_lp_solver requires at least A, b, and c as inputs.');
    end
    if size(A,1) ~= size(b,1)
        error('Dimension mismatch: A has %d rows, but b has %d elements.', size(A,1), size(b,1));
    end
    if size(A,2) ~= size(c,1)
        error('Dimension mismatch: A has %d columns, but c has %d elements.', size(A,2), size(c,1));
    end

    % Prepare MOSEK problem structure
    prob.c = c(:);                     % Ensure c is a column vector
    prob.a = sparse(A);                % Constraints matrix
    prob.blc = -inf*ones(size(b));     % No lower bounds on inequalities
    prob.buc = b;                      % Upper bounds of the inequalities
    prob.blx = [];                     % No lower bounds on x
    prob.bux = [];                     % No upper bounds on x

    % Default MOSEK options
    defaultOptions = struct('echo', 0);

    % Merge user-provided options
    if nargin >= 4 && isstruct(options)
        options = struct_merge(defaultOptions, options);
    else
        options = defaultOptions;
    end

    % Construct MOSEK command string
    cmd = sprintf('minimize echo(%d)', options.echo);

    % Solve the linear program using MOSEK
    [r, res] = mosekopt(cmd, prob);

    % Extract solution (handle cases where the solution might not exist)
    if isfield(res, 'sol') && isfield(res.sol, 'itr') && isfield(res.sol.itr, 'xx')
        x = res.sol.itr.xx(:);
    else
        x = NaN(size(c));  % Return NaNs if no solution is found
    end

    % Set flag similar to linprog
    flag = 0;
    if isfield(res.sol.itr, 'prosta')
        prosta = res.sol.itr.prosta;
        solsta = res.sol.itr.solsta;
        if strcmp(prosta, 'PRIMAL_AND_DUAL_FEASIBLE')
            flag = 1;
        elseif strcmp(prosta, 'DUAL_INFEASIBLE') && strcmp(solsta, 'DUAL_INFEAS')
            flag = -3;
        elseif strcmp(prosta, 'PRIMAL_INFEASIBLE') && strcmp(solsta, 'PRIMAL_INFEAS')
            flag = -3;  % Explicitly handle primal infeasibility
        elseif strcmp(solsta, 'UNKNOWN')
            warning('MOSEK returned an unknown solution status.');
        end
    else
        warning('Unexpected solver output: MOSEK did not return a proper status.');
    end 
end

function [x, flag] = scs_lp_solver(A, b, c, options)
    % SCS_LP_SOLVER Solves a linear program using the SCS solver
    %   Minimize:   c' * x
    %   Subject to: A * x = b
    %
    % Inputs:
    %   A        - Constraint matrix (m x n)
    %   b        - Right-hand side vector (m x 1)
    %   c        - Objective function coefficients (n x 1)
    %   options  - Struct with optional solver settings (optional)
    %
    % Outputs:
    %   x        - Solution vector
    %   flag     - Solver status: 
    %               1  -> Solved
    %              -3  -> Infeasible
    %              -2  -> Unbounded
    %               0  -> Other cases (failure)

    % Validate inputs
    if nargin < 3
        error('scs_lp_solver requires at least A, b, and c as inputs.');
    end
    if size(A,1) ~= size(b,1)
        error('Dimension mismatch: A has %d rows, but b has %d elements.', size(A,1), size(b,1));
    end
    if size(A,2) ~= size(c,1)
        error('Dimension mismatch: A has %d columns, but c has %d elements.', size(A,2), size(c,1));
    end

    % Ensure dense vectors if necessary
    if issparse(b), b = full(b); end
    if issparse(c), c = full(c); end

    % Prepare SCS input data
    data.A = A;     % Slack variables added
    data.c = c;     % Extend objective function
    data.b = b;

    % Define the positive cone
    K.l = size(A,1);

    % Default solver settings
    defaultSettings = struct('eps', 1e-8, 'verbose', 0);
    
    % Merge user-provided settings
    if nargin >= 4 && isstruct(options)
        settings = struct_merge(defaultSettings, options);
    else
        settings = defaultSettings;
    end
    
   % Solve LP using SCS
    [x, ~, ~, info] = scs_direct(data, K, settings); 

    % Map SCS status to flag values similar to linprog
    flag = 0;
    if strcmp(info.status, 'solved')
        flag = 1;
    elseif strcmp(info.status, 'infeasible')
        flag = -3;
    elseif strcmp(info.status, 'unbounded')
        flag = -2;
    else
        warning('Unexpected solver status in Newton reduction: %s', info.status);
    end
end


function [x, flag] = sedumi_lp_solver(A, b, c, options)
    % SEDUMI_LP_SOLVER Solves a linear program using the SeDuMi solver.
    %   Minimize:   c' * x
    %   Subject to: A * x <= b
    %
    % Inputs:
    %   A        - Constraint matrix (m x n)
    %   b        - Right-hand side vector (m x 1)
    %   c        - Objective function coefficients (n x 1)
    %   options  - (Optional) Struct with solver options
    %
    % Outputs:
    %   x        - Solution vector
    %   flag     - Solver status: 
    %               1  -> Solved
    %              -3  -> Unbounded
    %               0  -> Other cases (failure)

    % Validate inputs
    if nargin < 3
        error('sedumi_lp_solver requires at least A, b, and c as inputs.');
    end
    if size(A,1) ~= size(b,1)
        error('Dimension mismatch: A has %d rows, but b has %d elements.', size(A,1), size(b,1));
    end
    if size(A,2) ~= size(c,1)
        error('Dimension mismatch: A has %d columns, but c has %d elements.', size(A,2), size(c,1));
    end
    
    % Ensure column vectors
    b = b(:);
    c = c(:);

    % Reformulate problem for SeDuMi
    A_sedumi = [A speye(length(b))];        % Slack variables added
    c_sedumi = [c; zeros(length(b),1)];     % Extend objective function
    b_sedumi = b;

    % SeDuMi cone struct
    K.f = size(A,2);    % Free variables
    K.l = size(A,1);    % Linear inequalities

    % Default solver options
    defaultOptions = struct('fid', 0, 'eps', 1e-8);

    % Merge user-provided options if available
    if nargin >= 4 && isstruct(options)
        options = struct_merge(defaultOptions, options);
    else
        options = defaultOptions;
    end

    % Solve LP using SeDuMi
    [x, ~, info] = sedumi(A_sedumi, b_sedumi, c_sedumi, K, options); 

    % Extract relevant part of x (ignoring slack variables)
    x = x(1:size(A,2));

    % Set flag similar to linprog
    flag = 0;
    tolerance = 1e-2; % Tolerance for feasibility check

    if isfield(info, 'feasratio') && info.feasratio >= 1 - tolerance
        flag = 1;  % Feasible solution
    elseif isfield(info, 'pinf') && info.pinf == 0 && isfield(info, 'dinf') && info.dinf ~= 0
        flag = -3; % Unbounded case
    elseif isfield(info, 'numerr') && info.numerr > 1
        warning('SeDuMi encountered numerical issues (numerr > 1). Solution may be inaccurate.');
    end
end

function [x, flag] = linprog_solver(A, b, c, options)
    % LINPROG_SOLVER Solves a linear program using MATLAB's linprog.
    %   Minimize:   c' * x
    %   Subject to: A * x <= b
    %
    % Inputs:
    %   A        - Constraint matrix (m x n)
    %   b        - Right-hand side vector (m x 1)
    %   c        - Objective function coefficients (n x 1)
    %   options  - (Optional) Struct with solver options for linprog
    %
    % Outputs:
    %   x        - Solution vector
    %   flag     - Solver status:
    %               1  -> Solved successfully
    %              -3  -> Unbounded
    %               0  -> Other cases (failure)

    % Validate inputs
    if nargin < 3
        error('linprog_solver requires at least A, b, and c as inputs.');
    end
    if size(A,1) ~= size(b,1)
        error('Dimension mismatch: A has %d rows, but b has %d elements.', size(A,1), size(b,1));
    end
    if size(A,2) ~= size(c,1)
        error('Dimension mismatch: A has %d columns, but c has %d elements.', size(A,2), size(c,1));
    end

    % Ensure column vectors
    b = b(:);
    c = c(:);

    % Default solver options if none are provided
    if nargin < 4 || isempty(options)
        options = optimoptions('linprog', 'Display', 'off');
    end

    % Solve LP using linprog
    [x, ~, exitflag] = linprog(c, A, b, [], [], [], [], options);

    % Map linprog's exitflag to a custom flag
    switch exitflag
        case 1
            flag = 1;  % Optimal solution found
        case -3
            flag = -3; % Problem is unbounded
        otherwise
            flag = 0;  % Other failure cases
            warning('linprog failed with exitflag = %d. Solution may not be valid.', exitflag);
    end

end

function merged = struct_merge(defaults, custom)
    % Merges two structures, giving priority to the custom settings
    fields = fieldnames(custom);
    merged = defaults;
    for i = 1:numel(fields)
        merged.(fields{i}) = custom.(fields{i});
    end
end

