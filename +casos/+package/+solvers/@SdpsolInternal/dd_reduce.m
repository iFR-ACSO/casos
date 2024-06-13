function [sdp,args,opts] = dd_reduce(obj, sdp, opts)

args.dd_lbg = [];
args.dd_ubg = [];

% Verify DD cones in the constraints and create slack DD variables
if isfield(opts.Kc, 'ddm') && ~isempty(opts.Kc.ddm)
% Get DD cone sizes 
n_ddm = opts.Kc.ddm;

% get current number of constraints
n_g = length(sdp.g);

% update opts
if ~isfield(opts.Kc, 'lin')
    opts.Kc.lin = 0;
end    
if ~isfield(opts.Kx, 'ddm')
    opts.Kx.ddm = [];
end

% loop to iterate over each DD cone in the constraints
for i=1:length(n_ddm)
    % locate set of DD constraints (start from the end)
    dd_g = sdp.g(end - n_ddm(end-i+1)^2+1:end);

    % create slack variables that need to be DD
    dd_s = casadi.SX.sym(strjoin({'sg_', num2str(i)}, ''), n_ddm(end-i+1)^2, 1);

    % rewrite constraints and add to linear constraints
    dd_g = dd_g - dd_s;

    % add new constraints to the place of the linear
    % constraints
    sdp.g = [sdp.g(1:opts.Kc.lin); dd_g; sdp.g(opts.Kc.lin+1:end)];

    % remove the last constraints
    sdp.g = sdp.g(1: end - n_ddm(end-i+1)^2);

    % update lower and upper bound on the linear constraints
    % lower bounds is always zero
    args.dd_lbg = [args.dd_lbg; zeros(n_ddm(i)^2, 1)];
    % upper bounds is always infinity
    args.dd_ubg = [args.dd_ubg; zeros(n_ddm(i)^2, 1)];

    % update decision variables
    sdp.x = [sdp.x; dd_s];
    
    % add linear cones to constraints
    opts.Kc.lin = opts.Kc.lin + n_ddm^2;
    % add DD cones to decision variables
    opts.Kx.ddm = [opts.Kx.ddm; n_ddm];

end
% remove DD cone from constraints
opts.Kc = rmfield(opts.Kc, 'ddm');

end

if isfield(opts.Kx, 'ddm') && ~isempty(opts.Kx.ddm)
% Get DD cone sizes
n_ddm = opts.Kx.ddm;

% update opts
if ~isfield(opts.Kx, 'lin')
    opts.Kx.lin = 0;
end

for i=1:length(n_ddm)
    % Extract the ddm elements from sdp.x
    ddm_x = sdp.x(size(sdp.x,1)- n_ddm(i)^2+1:end);
    ddm_x = ddm_x.reshape(n_ddm(i), n_ddm(i)); % matrix of the form (n_ddm x n_ddm)

    % obtain diagonal elements of ddm_x
    ddm_x_diag = diag(ddm_x);
    
    % create slack variables (only corresponding to the upper
    % triangle except the diagonal elements
    if n_ddm(i)*(n_ddm(i)+1)/2 - n_ddm(i) >= 1
    ddm_s = casadi.SX.sym(strjoin({'sx_', num2str(i)}, ''), n_ddm(i)*(n_ddm(i)+1)/2 - n_ddm(i), 1);
    
    % create symmetric matrix with zero diagonal and upper and
    % bottom triangle equal to ddm_s            
    ind_s = triu(ones(n_ddm(i), n_ddm(i)),1); % extracts the upper triangle without the diagonal
    ind_d = ind_s + ind_s';
    
    % by the order of ddm_s, subs the 1 in ind_s with ddm_s
    s_matrix = casadi.SX(n_ddm(i)^2, 1);

    % find indexes where ind_s is 1
    s_matrix(find(ind_s(:))) = ddm_s;
    s_matrix = s_matrix.reshape(n_ddm(i), n_ddm(i));
    s_matrix = s_matrix + s_matrix';

    % add constraints of the type s_i - x_i >= 0 and s_i + x_i >= 0
    dd_g1 = s_matrix - ind_d.*ddm_x;
    dd_g2 = s_matrix + ind_d.*ddm_x;
    
    % create new constrains according to Gershgorin's circle theorem.
    dd_g3 = ddm_x_diag - sum(s_matrix,1)';
    
    % create constraints for symmetry of ddm_x
    dd_g4 = ddm_x - ddm_x';
    dd_g5 = ddm_x' - ddm_x;

    % new sdp.g
    n_new_g = size(sdp.g,1);
    sdp.g = [sdp.g(:); dd_g1(:); dd_g2(:); dd_g3(:); dd_g4(:); dd_g5(:)];
    n_new_g = size(sdp.g,1) - n_new_g;
    else
        ddm_s = [];
        n_new_g = 0;
    end

    % new sdp.x
    sdp.x = [sdp.x; ddm_s];

    
    
    % define size of new linear cones
    opts.Kx.lin = opts.Kx.lin + n_ddm(i)^2 + size(ddm_s,1); 
    opts.Kc.lin = opts.Kc.lin + n_new_g;

    % add new bounds
    % lower bounds is always zero
    args.dd_lbg = [args.dd_lbg; zeros(n_new_g,1)];
    % upper bounds is always infinity
    args.dd_ubg = [args.dd_ubg; inf(n_new_g,1)];
end

% remove ddm field from decision variables
opts.Kx = rmfield(opts.Kx, 'ddm');

end
