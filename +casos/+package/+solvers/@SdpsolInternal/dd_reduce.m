function [sdp,args,opts] = dd_reduce(obj, sdp, opts)
% Reduce a DD cone program to a linear program.

check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'dd');
check_cone(obj.get_cones,opts.Kc,'dd');

% initialize 
args.dd_lbg = [];
args.dd_ubg = [];

% get dimensions of linear and DD cones
Nl = get_dimension(obj.get_cones,opts.Kx,'lin');
Ml = get_dimension(obj.get_cones,opts.Kc,'lin');
Nd = get_dimension(obj.get_cones,opts.Kx,'dd');
Md = get_dimension(obj.get_cones,opts.Kc,'dd');

% Verify DD cones in the constraints and create slack DD variables
if isfield(opts.Kc,'dd')
    % loop to iterate over each DD cone in the constraints
    for i=1:length(Md)
        % locate set of DD constraints (start from the end)
        dd_g = sdp.g(end - Md(end-i+1)^2+1:end);
    
        % create slack variables that need to be DD
        dd_s = casadi.SX.sym(strjoin({'sg_', num2str(i)}, ''), Md(end-i+1)^2, 1);
    
        % rewrite constraints and add to linear constraints
        dd_g = dd_g - dd_s;
    
        % add new constraints to the place of the linear
        % constraints
        sdp.g = [sdp.g(1:opts.Kc.lin); dd_g; sdp.g(opts.Kc.lin+1:end)];
    
        % remove the last constraints
        sdp.g = sdp.g(1: end - Md(end-i+1)^2);
    
        % update lower and upper bound on the linear constraints
        % lower bounds is always zero
        args.dd_lbg = [args.dd_lbg; zeros(Md(i)^2, 1)];
        % upper bounds is always infinity
        args.dd_ubg = [args.dd_ubg; zeros(Md(i)^2, 1)];
    
        % update decision variables
        sdp.x = [sdp.x; dd_s];
    end
        
    % add linear cones to constraints
    Ml = Ml + sum(Md.^2);
    % add DD cones to decision variables
    Nd = [reshape(Nd,1,[]) reshape(Md,1,[])];

    % remove DD cone from constraints
    opts.Kc = rmfield(opts.Kc, 'dd');
end

if isfield(opts.Kx,'dd')
    % loop over each DD cone in the variables
    for i=1:length(Nd)
        % Extract the ddm elements from sdp.x
        dd_x = sdp.x(size(sdp.x,1)- Nd(i)^2+1:end);
        dd_x = dd_x.reshape(Nd(i), Nd(i)); % matrix of the form (n_ddm x n_ddm)
    
        % obtain diagonal elements of ddm_x
        ddm_x_diag = diag(dd_x);

        % number of elements in upper triangle minus the diagonal
        n_band = Nd(i)*(Nd(i)+1)/2 - Nd(i);
        
        if n_band >= 1
            % create slack variables 
            dd_s = casadi.SX.sym(strjoin({'sx_', num2str(i)}, ''), n_band, 1);
            
            % create symmetric matrix with zero diagonal and upper and
            % bottom triangle equal to ddm_s            
            ind_s = triu(ones(Nd(i), Nd(i)),1); % extracts the upper triangle without the diagonal
            ind_d = ind_s + ind_s';
            
            % find indexes where ind_s is 1
            S = casadi.Sparsity.nonzeros(Nd(i), Nd(i), find(ind_s));

            % by the order of ddm_s, subs the 1 in ind_s with ddm_s
            s_matrix = sparsity_cast(dd_s,S);
            % symmetrize
            s_matrix = s_matrix + s_matrix';
        
            % add constraints of the type s_i - x_i >= 0 and s_i + x_i >= 0
            dd_g1 = s_matrix - ind_d.*dd_x;
            dd_g2 = s_matrix + ind_d.*dd_x;
            
            % create new constrains according to Gershgorin's circle theorem.
            dd_g3 = ddm_x_diag - sum(s_matrix,1)';
            
            % create constraints for symmetry of ddm_x
            dd_g4 = dd_x - dd_x';
            dd_g5 = dd_x' - dd_x;
        
            % new sdp.g
            n_new_g = size(sdp.g,1);
            sdp.g = [sdp.g(:); dd_g1(:); dd_g2(:); dd_g3(:); dd_g4(:); dd_g5(:)];
            n_new_g = size(sdp.g,1) - n_new_g;

        else
            % nothing to do
            dd_s = [];
            n_new_g = 0;
        end
    
        % new sdp.x
        sdp.x = [sdp.x; dd_s];
        
        % define size of new linear cones
        Nl = Nl + Nd(i)^2 + size(dd_s,1); 
        Ml = Ml + n_new_g;
    
        % add new bounds
        % lower bounds is always zero
        args.dd_lbg = [args.dd_lbg; zeros(n_new_g,1)];
        % upper bounds is always infinity
        args.dd_ubg = [args.dd_ubg; inf(n_new_g,1)];
    end

    % remove dd field from decision variables
    opts.Kx = rmfield(opts.Kx, 'dd');
end

% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin',Nl);
opts.Kc = setfield(opts.Kc,'lin',Ml);

end
