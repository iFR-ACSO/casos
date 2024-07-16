function [sdp,args,map,opts] = sdd_reduce(obj, sdp, opts)
% Reduce a SDD cone program to SOCP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'sdd');
check_cone(obj.get_cones,opts.Kc,'sdd');

% initialize 
args.dd_lbg = [];
args.dd_ubg = [];

% get dimensions of linear and SDD cones
Nl = get_dimension(obj.get_cones,opts.Kx,'lin');
Ml = get_dimension(obj.get_cones,opts.Kc,'lin');
Nd = get_dimension(obj.get_cones,opts.Kx,'sdd');
Md = get_dimension(obj.get_cones,opts.Kc,'sdd');

% temporarily save sdp.x (original)
x_original = sdp.x;
g_original = sdp.g;


% transfer constraint on the x into the g
if Nd>0
    sdp.g = [sdp.g; sdp.x(end-Nd^2:end)];
    Md = Nd;
    Nd = 0;
    Nl = Nl + Md^2;
    opts.Kx = rmfield(opts.Kx, 'sdd');
end

% Verify SDD cones in the constraints and create slack SDD variables
if Md >0
    % loop to iterate over each SDD cone in the constraints
    for i=1:length(Md)

        % locate set of DD constraints (start from the end)
        dd_g = sdp.g(end - Md(end-i+1)^2+1:end);
        
        % matrix dimension
        n = Md(i);

        % Total number of pairs
        numPairs = (n-1)*n/2;

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

        % create casadi.SX.sym for the M
        M = casadi.SX.sym('M', 3*numPairs, 1);

        % obtain vec(Y)
        vecY = kronY*kron(speye(numPairs),M_selector); %*M;

        % method by equality constraints and casadi
        % dd_g-vecY*M == 0
        eq_constraints = dd_g-vecY*M;

        % constraints on M to be PSD
        ineq_constraints = kron(speye(numPairs),M_selector)*M; % all three in order shouldbelong to the psd
        
        % ToDo: 
        % add new constraints to the place of the linear
        % constraints
        sdp.g = [eq_constraints; ineq_constraints];
        %sdp.g = [sdp.g(1:opts.Kc.lin); dd_g; sdp.g(opts.Kc.lin+1:end)];
    
        % remove the last constraints
        % sdp.g = sdp.g(1: end - Md(end-i+1)^2);
    
        % update lower and upper bound on the linear constraints
        % lower bounds is always zero
        args.dd_lbg = [args.dd_lbg; zeros(length(eq_constraints), 1)];
        % upper bounds is always infinity
        args.dd_ubg = [args.dd_ubg; zeros(length(eq_constraints), 1)];
    
        % update decision variables
        sdp.x = [sdp.x; M];
    end
        
    % add linear cones to constraints
    Ml = Ml + length(sdp.g);
    Nl = Nl + length(M);

    % remove DD cone from constraints
    opts.Kc = rmfield(opts.Kc, 'sdd');
end

% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin',Nl);
opts.Kc = setfield(opts.Kc,'lin',length(eq_constraints));
opts.Kc = setfield(opts.Kc,'psd',repmat(2, 1, numPairs));

% map from new sdp.x to old
map.x = jacobian(x_original, sdp.x);
map.g = [speye(length(g_original)), sparse(length(g_original), length(sdp.g)-length(g_original))];

end

