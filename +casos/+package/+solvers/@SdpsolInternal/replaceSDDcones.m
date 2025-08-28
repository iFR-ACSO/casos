function [sdp, args, M_out, num_nlin, sdd_index, opts] = replaceSDDcones(~, sdp, sizes, Mlin, args, opts, field)
% Replace SDD cones (constraint or decision variable form)
%
% Inputs:
%   obj    - problem object with cone utilities
%   sdp    - struct containing decision variables (x) or constraints (g)
%   sizes  - vector of SDD cone sizes 
%            (Nsdd if field = 'x', Msdd if field = 'g')
%   Mlin   - number of existing linear constraints 
%            (relevant only if field = 'g'; ignored otherwise)
%   args   - struct containing bounds:
%              .sdd_lbg, .sdd_ubg  (for constraints)
%              .sdd_lbx, .sdd_ubx  (for variables)
%   opts   - struct with cone specifications (opts.Kc.sdd or opts.Kx.sdd)
%   field  - string flag:
%              'g' → replace SDD cones in constraints
%              'x' → replace SDD cones in decision variables
%
% Outputs:
%   sdp     - updated struct with SDD cones replaced by linear constraints
%   args    - updated struct with augmented bounds
%   M_out   - cell array of slack variable blocks (one per SDD cone)
%   num_nlin - number of equality and inequality constraints introduced
%   sdd_index - struct with indexes for the equality and inequality
%              constraints introduced
%   opts    - updated cone specification struct


nCones     = length(sizes);              % precompute number of SDD cones
total_eq   = nCones;                     % one equality per cone
total_ineq = sum((sizes.^2 - sizes)/2);  % total inequalities

% preallocate
M_out = cell(nCones,1);                  % slack variable blocks
eq_constraints   = cell(total_eq,1);     % equality constraints
ineq_constraints = cell(total_ineq,1);   % inequality constraints

% selector for PSD constraints on M (used later)
M_selector = sparse([1,2,3,4],[1,3,3,2],[1,1,1,1],4,3);

% counter for inequality constraints
ineq_counter = 1;

% loop over SDD cones
for i = 1:nCones
    idx = nCones - i + 1;   % process cones in reverse order
    n   = sizes(idx);       % size of current cone
    numPairs = (n-1)*n/2;   % number of off-diagonal pairs

    % locate block in sdp.(field)
    ind = [ length(sdp.(field)) + 1 - sum(sizes(end-i+1:end).^2);
            length(sdp.(field))     - sum(sizes(end-i+2:end).^2) ];
    sdd_block = sdp.(field)(ind(1):ind(2));

    % map M -> vec(Y)
    vecY = map_M_to_Y(n, numPairs);

    % slack variable
    M_out{i} = casadi.SX.sym([field '_M' num2str(i)], 3*numPairs, 1);

    % equality constraint
    eq_constraints{i} = sdd_block - vecY*M_out{i};

    % inequality constraints (PSD)
    for j = 1:numPairs
        ineq_constraints{ineq_counter} = M_selector*M_out{i}((j-1)*3+1:j*3);
        ineq_counter = ineq_counter + 1;
    end
end

% concatenate
eq_constraints   = vertcat(eq_constraints{:});      % combine equalities
ineq_constraints = vertcat(ineq_constraints{:});    % combine inequalities

% counts
num_eq   = length(eq_constraints);      % number of equality constraints
num_ineq = length(ineq_constraints);    % number of inequality constraints

% mapping indices
if Mlin ~= 0
    % constraints start after existing linear constraints
    start_idx = Mlin + 1;     
else
    start_idx = 1;
end
% indices of equality constraints
eq_idx   = start_idx : start_idx+num_eq-1;      
% indices of inequalities
ineq_idx = start_idx+num_eq : start_idx+num_eq+num_ineq-1; 

% insert constraints
if strcmp(field,'g')
    % remove old SDD
    sdp.g = sdp.g(1:end-sum(sizes.^2));  

    % insert new equality and inequality constraints
    if Mlin ~= 0
        sdp.g = [sdp.g(1:Mlin); eq_constraints; ineq_constraints; sdp.g(Mlin+1:end)];
    else
        sdp.g = [eq_constraints; ineq_constraints; sdp.g];
    end
else
    % for decision variables, same insertion but field='x'
    if Mlin ~= 0
        sdp.g = [sdp.g(1:Mlin); eq_constraints; ineq_constraints; sdp.g(Mlin+1:end)];
    else
        sdp.g = [eq_constraints; ineq_constraints; sdp.g];
    end
end

% save mapping
sdd_index.eq_idx   = eq_idx;
sdd_index.ineq_idx = ineq_idx;
sdd_index.num_eq   = num_eq;
sdd_index.num_ineq = num_ineq;  

% update bounds
args.dd_lbg = [args.dd_lbg; zeros(num_eq,1)];
args.dd_ubg = [args.dd_ubg; zeros(num_eq,1)];

num_nlin = num_eq + num_ineq;

% clean opts
if strcmp(field,'g') && isfield(opts.Kc,'sdd')
    opts.Kc = rmfield(opts.Kc,'sdd');
elseif strcmp(field,'x') && isfield(opts.Kx,'sdd')
    opts.Kx = rmfield(opts.Kx,'sdd');
end

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
