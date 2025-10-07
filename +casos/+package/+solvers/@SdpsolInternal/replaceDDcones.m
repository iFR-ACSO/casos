function [sdp, args, M_out, num_nlin, dd_index, opts] = replaceDDcones(~, sdp, sizes, Mlin, args, opts, field)
% Replace DD cones (constraint or decision variable form)
%
% Inputs:
%   obj    - problem object with cone utilities
%   sdp    - struct containing decision variables (x) or constraints (g)
%   sizes  - vector of DD cone sizes 
%            (Ndd if field = 'x', Mdd if field = 'g')
%   Mlin   - number of existing linear constraints 
%            (relevant only if field = 'g'; ignored otherwise)
%   args   - struct containing bounds:
%              .dd_lbg, .dd_ubg  (for constraints)
%              .dd_lbx, .dd_ubx  (for variables)
%   opts   - struct with cone specifications (opts.Kc.dd or opts.Kx.dd)
%   field  - string flag:
%              'g' → replace DD cones in constraints
%              'x' → replace DD cones in decision variables
%
% Outputs:
%   sdp      - updated struct with DD cones replaced by linear constraints
%   args     - updated struct with augmented bounds
%   M_out    - cell array of slack variable blocks (one per DD cone)
%   num_nlin - number of equality and inequality constraints introduced
%   dd_index - struct with indexes for the equality and inequality
%              constraints introduced
%   opts     - updated cone specification struct

    % add linear constraint to guarantee DD (used later)
    Mconstr_selector = sparse([1, 2, 3, 4, 1, 2, 3, 4],   ...
                              [1, 1, 2, 2, 3, 3, 3, 3],   ...
                              [1, 1, 1, 1, -1, 1, -1, 1], ...
                               4, 3);

    nCones = length(sizes);                  % precompute number of DD cones
    total_eq   = nCones;                     % one equality per cone
    total_ineq = sum((sizes.^2 - sizes)/2);  % total inequalities

    % preallocate
    M_out = cell(length(sizes),1);           % slack variable blocks
    eq_constraints   = cell(total_eq,1);     % equality constraints
    ineq_constraints = cell(total_ineq,1);   % inequality constraints

    % counter for inequality constraints
    ineq_counter = 1;  

    % loop to iterate over each DD cone (reverse order)
    for i = 1:nCones
        idx = nCones - i + 1;   % process cones in reverse order
        n   = sizes(idx);       % size of current cone
        numPairs = (n-1)*n/2;   % number of off-diagonal pairs

        % locate block in sdp.(field)
        ind = [ length(sdp.(field)) + 1 - sum(sizes(end-i+1:end).^2); ...
                length(sdp.(field))     - sum(sizes(end-i+2:end).^2) ];
        dd_block = sdp.(field)(ind(1):ind(2));

        % mapping
        vecY = map_M_to_Y(n, numPairs);

        % slack variable
        M_out{i} = casadi.SX.sym([field 'dd_M' num2str(i)], 3*numPairs, 1);

        % equality and inequality constraints
        eq_constraints{i} = dd_block - vecY*M_out{i};
        
        % inequality constraints
        for j = 1:numPairs
            ineq_constraints{ineq_counter} = Mconstr_selector * M_out{i}((j-1)*3+1:j*3);
            ineq_counter = ineq_counter + 1;
        end
    end

    % concatenate constraints
    eq_constraints   = vertcat(eq_constraints{:});   % combine equalities
    ineq_constraints = vertcat(ineq_constraints{:}); % combine inequalities

    % return counts
    num_eq   = length(eq_constraints);      % number of equality constraints
    num_ineq = length(ineq_constraints);    % number of inequality constraints

    % build mapping indices
    if Mlin ~= 0
        % position where eq_constraints start
        start_idx = Mlin + 1;                    
        eq_idx    = start_idx : start_idx+num_eq-1;
        ineq_idx  = start_idx+num_eq : start_idx+num_eq+num_ineq-1;
    else
        start_idx = 1;
        eq_idx    = start_idx : start_idx+num_eq-1;
        ineq_idx  = start_idx+num_eq : start_idx+num_eq+num_ineq-1;
    end

    % Special handling depending on field
    if strcmp(field,'g')
        % remove previous constraints related to dd
        sdp.g = sdp.g(1:end-sum(sizes.^2));

        % insert new constraints
        if Mlin ~= 0
            sdp.g = [sdp.g(1:Mlin); eq_constraints; ineq_constraints; sdp.g(Mlin+1:end)];
        else
            sdp.g = [eq_constraints; ineq_constraints; sdp.g];
        end
    else
        % for decision variables, just add constraints to sdp.g
        if Mlin ~= 0
            sdp.g = [sdp.g(1:Mlin); eq_constraints; ineq_constraints; sdp.g(Mlin+1:end)];
        else
            sdp.g = [eq_constraints; ineq_constraints; sdp.g];
        end
    end

    % save mapping for later retrieval (duals)
    dd_index.eq_idx   = eq_idx;    % indices in sdp.g for eq_constraints
    dd_index.ineq_idx = ineq_idx;  % indices in sdp.g for ineq_constraints
    dd_index.num_eq   = num_eq;
    dd_index.num_ineq = num_ineq;

    % update bounds
    args.dd_lbg = [args.dd_lbg; zeros(num_eq,1); zeros(num_ineq,1)];
    args.dd_ubg = [args.dd_ubg; zeros(num_eq,1);  inf(num_ineq,1)];

    % retrieve the number of new constraints
    num_nlin = num_eq + num_ineq;

    % clean up opts
    if strcmp(field,'g') && isfield(opts.Kc,'dd')
        opts.Kc = rmfield(opts.Kc,'dd');
    elseif strcmp(field,'x') && isfield(opts.Kx,'dd')
        opts.Kx = rmfield(opts.Kx,'dd');
    end
end


function vecY = map_M_to_Y(n, numPairs)
    % build vectorized Y
    M_selector = sparse([1, 4, 2, 3], [1, 2, 3, 3], [1, 1, 1, 1], 4, 3);

    % build indexes for T and S
    i_pair = ceil((2*n - 1 - sqrt((2*n - 1)^2 - 8*(1:numPairs))) / 2);
    j_pair = (1:numPairs) - (i_pair-1).*n + i_pair.*(i_pair-1)./ 2 + i_pair;

    % avoid using this for loop
    T = arrayfun(@(k) sparse([1, 2], [i_pair(k), j_pair(k)], [1, 1], 2, n), 1:numPairs, 'UniformOutput', false);
    S = arrayfun(@(k) sparse([i_pair(k), j_pair(k)], [1, 2], [1, 1], n, 2), 1:numPairs, 'UniformOutput', false);

    % selects from [M1 M3; M3 M2] where to put in matrix vec(Y)
    kronY = cell2mat(arrayfun(@(k) kron(T{k}', S{k}), 1:numPairs, 'UniformOutput', false));

    % obtain vec(Y)
    vecY = kronY*kron(speye(numPairs),M_selector); %*M;
end
