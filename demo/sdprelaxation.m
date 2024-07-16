% Test and compare SDP with DD and SDD relaxation

% set seed
rng(0)

% dimension of SDP
ndim = 5;

% decision variables
x = casadi.SX.sym('x', 2, 1);

% create random SDP data
A = randn(ndim, ndim);
B = randn(ndim, ndim);
A = A+A';
B = B+B';

% define SDP problem
sdp.f = x(1);
sdp.g = eye(ndim) + x(1)*A + x(2)*B;
sdp.g = sdp.g(:);
sdp.x = x;

opts = struct;

%% Solve SDP with psd constraint
% define cones
opts.Kx = struct('lin', 2);
opts.Kc = struct('psd', ndim);

% initialize solver
S = casos.sdpsol('S','mosek',sdp,opts);

% solve with equality constraints
tic
sol = S('lbx',-inf,'ubx',+inf);
elapsed_time = toc;
if S.stats.UNIFIED_RETURN_STATUS == "SOLVER_RET_SUCCESS"
    fprintf('Elasped time with psd: %d \n',elapsed_time);
    fprintf('Solution with psd: %d \n', sol.f.full);
    sol_psd_x = full(sol.x);
else
    fprintf('Feasibility issues while solving with psd');
end

%% Solve SDP with relaxation to DD 

% define cones
opts.Kx = struct('lin', 2);
opts.Kc = struct('dd', ndim);

% initialize solver
S = casos.sdpsol('S','mosek',sdp,opts);

% solve with equality constraints
tic
sol = S('lbx',-inf,'ubx',+inf);
elapsed_time = toc;
if S.stats.UNIFIED_RETURN_STATUS == "SOLVER_RET_SUCCESS"
    fprintf('Elapsed time with dd: %d \n', elapsed_time);
    fprintf('Solution with dd: %d \n', sol.f.full);
    sol_dd_x = full(sol.x);
else
    fprintf('Feasibility issues while solving with dd');
end

%% Solve SDP with relaxation to SDD

% define cones
opts.Kx = struct('lin', 2);
opts.Kc = struct('sdd', ndim);

% initialize solver
S = casos.sdpsol('S','mosek',sdp,opts);

% solve with equality constraints
tic
sol = S('lbx',-inf,'ubx',+inf);
elapsed_time = toc;
if S.stats.UNIFIED_RETURN_STATUS == "SOLVER_RET_SUCCESS"
    fprintf('Elapsed time with sdd: %d \n', elapsed_time);
    fprintf('Solution with sdd: %d \n', sol.f.full);
    sol_sdd_x = full(sol.x);
else
    fprintf('Feasibility issues while solving with sdd');
end

%% plot the set of feasible solutions for the SDP with PSD, DD and SDD

% create grid 
xp = -0.5:0.002:0.5;
yp = -0.5:0.002:0.5;
[xp, yp] = meshgrid(xp,yp);
xp = xp(:);
yp = yp(:);

% create variables to store feasible points
fplot = cell(3);

for k = 1:length(xp)
    
    % verify if it is psd
    if all(eig(eye(ndim) + xp(k)*A + yp(k)*B)>=0)
        fplot{1} = [fplot{1}; xp(k) yp(k)];
    end
    
    % verify if it is diagonally dominant
    if isDD(eye(ndim) + xp(k)*A + yp(k)*B)
        fplot{2} = [fplot{2}; xp(k) yp(k)];
    end

    % verify if it is scaled diagonally dominant
    if isSDD(eye(ndim) + xp(k)*A + yp(k)*B)
        fplot{3} = [fplot{3}; xp(k) yp(k)];
    end

end

% initialize figure
figure(1)
legend('show');
xlabel('x_1'); ylabel('x_2')
hold on

legend_names = {'PSD feasible set', 'DD feasible set', 'SDD feasible set'};
% plot feasible set with psd, dd and sdd
for i = 1:3
    [k,~] = convhull(fplot{i});
    plot(fplot{i}(k,1),fplot{i}(k,2), 'DisplayName', legend_names{i})
end

% if solutions exist, plot them
legend_names = {'PSD optimal solution', 'DD optimal solution', 'SDD optimal solution'};
sol_names = {'sol_psd_x', 'sol_dd_x', 'sol_sdd_x'};

for i=1:3
    if exist(sol_names{i}, 'var')
        aux = eval(sol_names{i});
        scatter(aux(1), aux(2), 'DisplayName', legend_names{i});
    end
end

%%
% Function that verifies if matrix is diagonally dominat 
function out = isDD(A)
    
    [m, n] = size(A);

    % Check if the matrix is square
    if m ~= n
        error('Matrix must be square');
    end

    if ~(all(diag(A) >= 0) && issymmetric(A))
        out = 0;
        return
    end
    diagA = diag(diag(A));
    out = all(sum(2*diagA-abs(A),2) >= 0);
end

% Function that verifies if matrix is scaled diagonally dominat 
function out = isSDD(A)
    
    % Get size of A
    [m, n] = size(A);

    % Check if the matrix is square
    if m ~= n
        error('Matrix must be square');
    end
    
    % Verify if matrix is symmetric and all diagonal elements are
    % nonnegative
    if ~(all(diag(A) >= 0) && issymmetric(A))
        out = 0;
        return
    end

    % Set up and solve the linear programming problem to find scaling vector d
    f = -ones(m, 1);    % Objective function: minimize -sum(d)
    Aineq = [];         % Inequality constraints
    bineq = [];         % Inequality constraints
    Aeq = [];           % Equality constraints
    beq = [];           % Equality constraints
    lb = zeros(m, 1);   % Lower bounds for d: d > 0
    ub = [];            % No upper bounds

    % Construct the inequality constraints for diagonal dominance
    for i = 1:m
        row = -abs(A(i, :) / A(i, i));  % Compute the entire row, including the diagonal
        row(i) = 1;                     % Set the diagonal element to 1
        Aineq = [Aineq; row];
        bineq = [bineq; 1]; 
    end

    % Solve the linear programming problem
    options = optimoptions('linprog', 'Display', 'none');
    [~, ~, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);
    out = (exitflag==1);
    
end


