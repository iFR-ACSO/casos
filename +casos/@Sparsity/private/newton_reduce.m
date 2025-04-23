function Lz = newton_reduce(Pdegmat,Zdegmat,solver)
% Removes monomials outside half Newton polytope
% Strategy inpired from: Simplification Methods for Sum-of-Squares 
% Programs, Peter Seiler et al.

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
        % for the conic solver:
        lin = struct('g', casadi.Sparsity.dense(size(c,1),size(c,2)), 'a', casadi.Sparsity.dense(size(A)));

        % Set options for solver
        opts.Kx = struct('lin', size(A,2));
        opts.Kc = struct('lin', length(b));
        opts.error_on_fail = false;

        % Build the linear program
        Slin = casos.conic('S', solver, lin, opts);
    end

    % Solve and obtain the solution
    sol = Slin('h',    sparse(size(A,2), size(A,2)), ...
                'g',    c,   ...
                'a',    A,   ...
                'lba',  -inf,    ...        
                'uba',  b, ...    
                'cba',  [], ...
                'lbx',  -inf, ...        
                'ubx',  +inf, ...
                'cbx',  [], ...
                'x0',   sparse(size(A,2),1), ...
                'lam_a0', sparse(length(b),1) ,...
                'lam_x0', sparse(length(c),1));

    % Extract the solution
    x = full(sol.x);

    switch Slin.stats.UNIFIED_RETURN_STATUS
        case 'SOLVER_RET_SUCCESS'
            flag = 1;
        case 'SOLVER_RET_LIMITED'
            flag = 0;
        case 'SOLVER_RET_INFEASIBLE'
            flag = -3;
            x = zeros(length(c),1);
        case 'SOLVER_RET_NAN'
            flag = -4;
        case 'SOLVER_RET_UNKNOWN'
            flag = -1;
        otherwise
            flag = -99;
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
