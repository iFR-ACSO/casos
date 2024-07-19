function [sdp,args,map,opts] = dd_reduce(obj, sdp, opts)
% Reduce a DD cone program to LP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'dd');
check_cone(obj.get_cones,opts.Kc,'dd');

% initialize 
% ToDo: due to ordering, create a map for this 
args.dd_lbg = [];
args.dd_ubg = [];

% get dimensions of cones in the program decision variables
Nlin = get_dimension(obj.get_cones,opts.Kx,'lin');
Nlor = get_dimension(obj.get_cones,opts.Kx,'lor');
Nrot = get_dimension(obj.get_cones,opts.Kx,'rot');
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd');
Ndd  = get_dimension(obj.get_cones,opts.Kx,'dd');
Nsdd = get_dimension(obj.get_cones,opts.Kx,'sdd');

% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin');
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor');
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot');
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd');
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');
Msdd = get_dimension(obj.get_cones,opts.Kc,'sdd');

% temporarily save sdp.x (original)
x_original = sdp.x;
g_original = sdp.g;

% add linear constraint to guarantee DD
Mconstr_selector = sparse([1, 2, 3, 4, 1, 2, 3, 4],   ...
                          [1, 1, 2, 2, 3, 3, 3, 3],   ...
                          [1, 1, 1, 1, -1, 1, -1, 1], ...
                          4, 3);

% create zero initial 
num_eq_g = 0;
num_eq_x = 0;
num_ineq_g = 0;
num_ineq_x = 0;

% Verify DD cones in the constraints and create slack DD variables
M_g = cell(length(Mdd),1);
if isfield(opts.Kc,'dd')
    % loop to iterate over each DD cone in the constraints
    eq_constraints_g   = cell(length(Mdd),1);
    ineq_constraints_g = cell(length(Mdd),1);

    for i=1:length(Mdd)
        % locate set of DD constraints (start from the end)
        ind = [length(sdp.g)+1-sum(Mdd(end-i+1:end).^2); length(sdp.g)-sum(Mdd(end-i+2:end).^2)]; 

        % locate set of DD constraints (start from the end)
        dd_g = sdp.g(ind(1):ind(2));
        
        % matrix dimension
        n = Mdd(end-i+1);

        % Total number of pairs
        numPairs = (n-1)*n/2;

        % create map from M to vec(Y)
        vecY = map_M_to_Y(n, numPairs);

        % create casadi.SX.sym for the M
        M_g{i} = casadi.SX.sym(['M_g' num2str(i)], 3*numPairs, 1);

        % equality constraints
        eq_constraints_g{i} = dd_g-vecY*M_g{i};  

        % inequality constraints
        ineq_constraints_g{i} = kron(speye(numPairs), Mconstr_selector)*M_g{i};
    end

    % remove previous constraints related to dd from sdp.g
    sdp.g = sdp.g(1:end-sum(Mdd.^2));

    % add new constraints to the place of the linear constraints
    sdp.g = [sdp.g; vertcat(eq_constraints_g{:}); vertcat(ineq_constraints_g{:})];

    num_eq_g   = length(vertcat(eq_constraints_g{:}));
    num_ineq_g = length(vertcat(ineq_constraints_g{:}));

    % update lower and upper bound on the linear constraints
    args.dd_lbg = [zeros(num_eq_g, 1); zeros(num_ineq_g, 1)];
    
    % upper bounds is always infinity
    args.dd_ubg = [zeros(num_eq_g, 1); inf(num_ineq_g,1)];

    % remove DD cone from constraints
    opts.Kc = rmfield(opts.Kc, 'dd');
end

% Verify DD cones in decision variables
M_x = cell(length(Ndd),1);
if isfield(opts.Kx, 'dd')
    
    % loop to iterate over each DD cone decision variables
    eq_constraints_x   = cell(length(Ndd),1);
    ineq_constraints_x = cell(length(Ndd),1);

    for i=1:length(Ndd)
        % locate set of DD constraints (start from the end)
        ind = [length(sdp.x)+1-sum(Ndd(end-i+1:end).^2); length(sdp.x)-sum(Ndd(end-i+2:end).^2)]; 

        % locate set of DD constraints (start from the end)
        dd_x = sdp.x(ind(1):ind(2));
        
        % matrix dimension
        n = Ndd(end-i+1);

        % Total number of pairs
        numPairs = (n-1)*n/2;

        % create map from M to vec(Y)
        vecY = map_M_to_Y(n, numPairs);

        % create casadi.SX.sym for the M
        M_x{i} = casadi.SX.sym(['M_x' num2str(i)], 3*numPairs, 1);

        % method by equality constraints and casadi
        eq_constraints_x{i} = dd_x-vecY*M_x{i};

        % save variables and constraints
        ineq_constraints_x{i} = kron(speye(numPairs), Mconstr_selector)*M_x{i};
    end

    % add new constraints to the place of the linear constraints
    sdp.g = [sdp.g; vertcat(eq_constraints_x{:}); vertcat(ineq_constraints_x{:})];

    % number of new equality and cone constainemnt constraints
    num_eq_x   = length(vertcat(eq_constraints_x{:}));
    num_ineq_x = length(vertcat(ineq_constraints_x{:}));

    % update lower and upper bound on the linear constraints
    args.dd_lbg = [args.dd_lbg; zeros(num_eq_x, 1); zeros(num_ineq_x, 1)];

    % upper bounds is always infinity
    args.dd_ubg = [args.dd_ubg; zeros(num_eq_x, 1); inf(num_ineq_x,1)];

    % remove DD cone from constraints
    opts.Kx = rmfield(opts.Kx, 'dd');
end

% update decision variables
sdp.x = [sdp.x; vertcat(M_g{:}); vertcat(M_x{:})];

% add linear cones to constraints
Mlin = Mlin + num_ineq_x + num_eq_x + num_eq_g + num_ineq_g;

% add the variables M to sdp.x
Nlin = Nlin + sum(Ndd.^2) + length(vertcat(M_g{:})) + length(vertcat(M_x{:}));

% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin',Nlin);
opts.Kc = setfield(opts.Kc,'lin',Mlin);

% map from new sdp.x to old
map.x = jacobian(x_original, sdp.x);
map.g = [speye(length(g_original)), sparse(length(g_original), length(sdp.g)-length(g_original))];

end



function vecY = map_M_to_Y(n, numPairs)
    % build vectorized Y
    M_selector = sparse([1, 4, 2, 3], [1, 2, 3, 3], [1, 1, 1, 1], 4, 3);

    % build indexes for T and S
    i_pair = ceil((2*n - 1 - sqrt((2*n - 1)^2 - 8*(1:numPairs))) / 2);
    j_pair = (1:numPairs) - (i_pair-1).*n + i_pair.*(i_pair-1)./ 2 + i_pair;

    % Avoid using this for loop
    T = arrayfun(@(k) sparse([1, 2], [i_pair(k), j_pair(k)], [1, 1], 2, n), 1:numPairs, 'UniformOutput', false);
    S = arrayfun(@(k) sparse([i_pair(k), j_pair(k)], [1, 2], [1, 1], n, 2), 1:numPairs, 'UniformOutput', false);

    % selects from [M1 M3; M3 M2] where to put in matrix vec(Y)
    kronY = cell2mat(arrayfun(@(k) kron(T{k}', S{k}), 1:numPairs, 'UniformOutput', false));

    % obtain vec(Y)
    vecY = kronY*kron(speye(numPairs),M_selector); %*M;
end





