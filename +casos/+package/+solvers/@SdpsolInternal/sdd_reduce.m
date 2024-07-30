function [sdp,args,map,opts] = sdd_reduce(obj, sdp, opts, args)
% Reduce a SDD cone program to SOCP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'sdd');
check_cone(obj.get_cones,opts.Kc,'sdd');

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

M_g = cell(length(Msdd),1);
% loop to iterate over each SDD cone in the constraints
eq_constraints_g   = cell(length(Msdd),1);
ineq_constraints_g = cell(length(Msdd),1);

% Verify SDD cones in the constraints and create slack SDD variables
if isfield(opts.Kc, 'sdd')

    for i=1:length(Msdd)
        % locate set of SDD constraints (start from the end)
        ind = [length(sdp.g)+1-sum(Msdd(end-i+1:end).^2); length(sdp.g)-sum(Msdd(end-i+2:end).^2)]; 

        % locate set of SDD constraints (start from the end)
        dd_g = sdp.g(ind(1):ind(2));
        
        % matrix dimension
        n = Msdd(end-i+1);

        % Total number of pairs
        numPairs = (n-1)*n/2;

        vecY = map_M_to_Y(n, numPairs);
        
        % create casadi.SX.sym for the M
        M_g{i} = casadi.SX.sym(['M_g' num2str(i)], 3*numPairs, 1);

        % method by equality constraints and casadi
        eq_constraints_g{i} = dd_g-vecY*M_g{i};

        M_selector = sparse([1, 2, 3, 4], [1, 3, 3, 2], [1, 1, 1, 1], 4, 3);

        % constraints on M to be PSD
        ineq_constraints_g{i} = kron(speye(numPairs),M_selector)*M_g{i};
    end

    % remove previous constraints related to dd from sdp.g
    sdp.g = sdp.g(1:end-sum(Msdd.^2));

    % add new constraints to the place of the linear constraints
    if Mlin~=0
        sdp.g = [sdp.g(1:Mlin); vertcat(eq_constraints_g{:}); sdp.g(Mlin+1:end)];
    else
        sdp.g = [vertcat(eq_constraints_g{:}); sdp.g];
    end
    sdp.g = [sdp.g; vertcat(ineq_constraints_g{:})];

    % get number of equality and inequality constraints
    num_eq   = length(vertcat(eq_constraints_g{:}));
    
    % update lower and upper bound on the linear constraints
    args.dd_lbg = [args.dd_lbg; zeros(num_eq, 1)];
    
    % upper bounds is always infinity
    args.dd_ubg = [args.dd_ubg; zeros(num_eq, 1)];

    % remove DD cone from constraints
    opts.Kc = rmfield(opts.Kc, 'sdd');
end

% Verify SDD cones in the decision variables
M_x = cell(length(Nsdd),1);
% loop to iterate over each SDD cone in the constraints
eq_constraints_x   = cell(length(Nsdd),1);
ineq_constraints_x = cell(length(Nsdd),1);
if isfield(opts.Kx, 'sdd')

    for i=1:length(Nsdd)
        % locate set of SDD constraints (start from the end)
        ind = [length(sdp.x)+1-sum(Nsdd(end-i+1:end).^2); length(sdp.x)-sum(Nsdd(end-i+2:end).^2)]; 

        % locate set of SDD constraints (start from the end)
        dd_x = sdp.x(ind(1):ind(2));
        
        % matrix dimension
        n = Nsdd(end-i+1);

        % Total number of pairs
        numPairs = (n-1)*n/2;

        % create map from M to vec(Y)
        vecY = map_M_to_Y(n, numPairs);

        % create casadi.SX.sym for the M
        M_x{i} = casadi.SX.sym(['M_x' num2str(i)], 3*numPairs, 1);

        % method by equality constraints and casadi
        eq_constraints_x{i} = dd_x-vecY*M_x{i};

        % build vectorized Y
        M_selector = sparse([1, 2, 3, 4], [1, 3, 3, 2], [1, 1, 1, 1], 4, 3);

        % constraints on M to be PSD
        ineq_constraints_x{i} = kron(speye(numPairs),M_selector)*M_x{i}; 
    end

    % add new constraints to the place of the linear constraints
    if Mlin~=0
        sdp.g = [sdp.g(1:Mlin); vertcat(eq_constraints_x{:}); sdp.g(Mlin+1:end)];
    else
        sdp.g = [vertcat(eq_constraints_x{:}); sdp.g];
    end
    sdp.g = [sdp.g; vertcat(ineq_constraints_x{:})];
    
    % number of new equality and cone constainemnt constraints
    num_eq   = length(vertcat(eq_constraints_x{:}));

    % update lower and upper bound on the linear constraints
    args.dd_lbg = [args.dd_lbg; zeros(num_eq, 1)];

    % upper bounds is always infinity
    args.dd_ubg = [args.dd_ubg; zeros(num_eq, 1)];

    % remove DD cone from constraints
    opts.Kx = rmfield(opts.Kx, 'sdd');
end

% update decision variables
sdp.x = [sdp.x(1:Nlin); vertcat(M_g{:}); vertcat(M_x{:}); sdp.x(Nlin+1:end)];

% update lbx and ubx 
args.dd_lbx = [args.dd_lbx; -inf(sum(Nsdd.^2) + length(vertcat(M_g{:})) + length(vertcat(M_x{:})),1)];
args.dd_ubx = [args.dd_ubx; inf(sum(Nsdd.^2) + length(vertcat(M_g{:})) + length(vertcat(M_x{:})),1)];


% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin', Nlin + sum(Nsdd.^2) + length(vertcat(M_g{:})) + length(vertcat(M_x{:})));
opts.Kc = setfield(opts.Kc,'lin', Mlin + length(vertcat(eq_constraints_x{:}))+length(vertcat(eq_constraints_g{:})));
opts.Kc = setfield(opts.Kc,'psd', [Mpsd, repmat(2, 1, length(vertcat(M_g{:}))/3+length(vertcat(M_x{:}))/3)]);

% map from new sdp.x to old
map.x = jacobian(x_original, sdp.x);
map.g = [speye(length(g_original)), sparse(length(g_original), length(sdp.g)-length(g_original))];

end


function vecY = map_M_to_Y(n, numPairs)
    % build vectorized Y
    M_selector = sparse([1, 2, 3, 4], [1, 3, 3, 2], [1, 1, 1, 1], 4, 3);

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





