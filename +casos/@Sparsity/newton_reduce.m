function Lz = newton_reduce(obj, Pdegmat, Zdegmat, solver)
% Removes monomials outside half Newton polytope
% Strategy inpired from: Simplification Methods for Sum-of-Squares 
% Programs, Peter Seiler et al.

% list of available solvers for the newton polytope reduction
list_solvers = {'linprog', 'sedumi', 'mosekopt'};
[solver_available, solver_id] = ismember(solver,list_solvers);

if solver_available==0
    error('Specified solver is not available.\n');
end

switch solver_id
    case 1  % options for linprog (Optimization toolbox from Matlab)
        options = optimoptions('linprog');
        options.Display = 'off';   
        lp_solver = @(A,b,c,options) linprog_solver(A,b,c,options);
    case 2  % options for Sedumi
        options.fid = 0;
        lp_solver = @(A,b,c,options) sedumi_lp_solver(A,b,c,options);
    case 3  % options for Mosek
        options = struct();
        lp_solver = @(A,b,c,options) mosek_lp_solver(A,b,c,options);
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
    prob.c = c;
    prob.a = sparse(A);             % Constraints matrix
    prob.blc = -inf*ones(size(b));  % No lower bounds on inequalities
    prob.buc = b;                   % Upper bounds of the inequalities
    prob.blx = [];                  % Lower bounds on x
    prob.bux = [];                  % No upper bounds

    % Solve the linear program using MOSEK
    [r, res] = mosekopt('minimize echo(0)', prob);

    x = res.sol.itr.xx(:);

    % set flag similar to linprog
    flag = 0;
    if res.sol.itr.prosta == "PRIMAL_AND_DUAL_FEASIBLE"
        flag = 1;
    end
    if res.sol.itr.prosta == "DUAL_INFEASIBLE" && res.sol.itr.solsta =="DUAL_INFEAS"
        flag =-3;
    end
    
end


function [x, flag] = sedumi_lp_solver(A, b, c, options)

    A_sedumi = [A speye(length(b))];
    c_sedumi = [c; zeros(length(b),1)];
    b_sedumi = b;
    K.f = size(A,2);
    K.l = size(A,1);
    options.fid = 0;
    [x, ~, info] = sedumi(A_sedumi, b_sedumi, c_sedumi, K, options); 
    x = x(1:size(A,2));

    % set flag similar to linprog
    flag = 0;
    if info.feasratio >= 1-0.01
        flag = 1;   % feasible solution
    end
    if info.pinf==0 && info.dinf~=0
        flag = -3;  % unbounded case
    end

end

function [x, flag] = linprog_solver(A,b,c,options)
    [x,~,flag] = linprog(c, A, b, [], [], [], [], options); 
end

