function Lz = newton_reduce(Pdegmat,Zdegmat,solver)
% Removes monomials outside half Newton polytope
% Strategy inpired from: Simplification Methods for Sum-of-Squares 
% Programs, Peter Seiler et al.

% list of available solvers for the newton polytope reduction
list_solvers = {'sedumi', 'mosek', 'scs'};
[solver_available, ~] = ismember(solver,list_solvers);

if solver_available==0
    error('Specified solver is not available.\n');
end

% build LP (part 1)
bfixed = [zeros(size(Pdegmat,1),1); 1];

% trivial removal of monomials
keep_trivial = ismember(Zdegmat*2,Pdegmat,'rows');
keep = true(size(keep_trivial));

% initialize some parameters for optimization
p_A = [];

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
    
    A = -F_struc(1:end,2:end);  % Linear constraint matrix
    b = F_struc(1:end,1);       % Right-hand side vector

    % Build LP if dimensions mismatch or at first iteration 
    % (reuse previous LP)
    if i==1 || any(size(p_A) ~=size(A))
        % Define symbolic variables for LP (parameters for LP)
        p_A = casadi.SX.sym('p_A', size(A));    % Symbolic for A
        p_c = casadi.SX.sym('p_c', size(c));    % Symbolic for cost vector c

        % Define decision variable vector
        newton_x = casadi.SX.sym('x', size(A,2), 1);

        % Define linear program components
        lin.f = p_c'*newton_x;          % Objective function
        lin.x = newton_x;               % Decision variables
        lin.g = p_A*newton_x;           % Constraints
        lin.p = [p_A(:); p_c(:)];       % Parameters (flattened)

        % Set options for solver
        opts.Kx = struct('lin', size(A,2));
        
        % Build the linear program
        S = casos.sdpsol('S', solver, lin, opts);
    end

    % Solve and obtain the solution
    sol = S('lbg', -inf(length(b),1), ...       % Lower bounds for constraints
            'ubg', b, ...                       % Upper bounds for constraints
            'lbx', -inf, ...                    % Lower bounds for decision variables
            'ubx', +inf, ...                    % Upper bounds for decision variables
            'p', [full(A(:)); full(c(:))]);     % Parameters (flattened)

    % Extract the solution
    x = full(sol.x);

    % Get the status code from solver
    flag = mapSolverReturn(S.stats.UNIFIED_RETURN_STATUS);

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


% Maps solver return values to numerical codes and descriptions.
function [code, description] = mapSolverReturn(solver_ret)
    
    mapping = struct( ...
        'SOLVER_RET_SUCCESS',    {1,  'Function converged to a solution.'}, ...
        'SOLVER_RET_LIMITED',    {0,  'Solution is feasible but does not fully satisfy tolerances.'}, ...
        'SOLVER_RET_INFEASIBLE', {-3, 'No feasible point was found or is unbounded'}, ...
        'SOLVER_RET_NAN',        {-4, 'NaN value encountered during execution.'}, ...
        'SOLVER_RET_UNKNOWN',    {-1, 'Solver did not return a clear result.'} ...
    );

    if isfield(mapping, solver_ret)
        values = {mapping.(solver_ret)};  % Extract the cell array
        code = values{1};                 % First element: numerical code
        description = values{2};          % Second element: description
    else
        code = -99;
        description = 'Unknown solver return value.';
    end
end