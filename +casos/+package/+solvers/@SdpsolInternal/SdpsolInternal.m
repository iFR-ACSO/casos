classdef SdpsolInternal < casos.package.solvers.SolverCallback & matlab.mixin.Copyable
% Internal interface for convex cone (SDP) solvers.

properties (Access=private)
    solver;
end

properties (Access=protected)
    fhan;
    ghan;
end

methods (Access=private)
    buildproblem(obj,prob,data,opts);
end

methods (Static)
    function cones = get_cones
        % Return supported cones.
        cones = casos.package.solvers.ConicSolver.get_cones;
    end
end

methods
    argout = eval(obj,argin);

    function obj = SdpsolInternal(name,solver,sdp,opts)
        obj@casos.package.solvers.SolverCallback;

        if isa(name,'casos.package.solvers.SdpsolInternal')
            % copy constructor
            obj.solver = name.solver;

            if nargin > 1
                % change input function
                obj.fhan = solver;
            else
                obj.fhan = name.fhan;
            end
            obj.ghan = name.ghan;

            % construct CasADi callback
            construct(obj,name.name);

            return
        end

        % options
        if nargin < 4
            opts = struct;
        else
            if isfield(opts, 'Kx') && isfield(opts.Kx, 'ddm')
                [sdp, opts] = dd_relax(obj, sdp, opts);
            end
        end

        % decision variables
        x = sdp.x;
        % parameter
        if isfield(sdp,'p')
            p = sdp.p;
        else
            p = [];
        end
        % constraint function (vectorized)
        sdp_g = sdp.g(:);

        % quadratic cost
        H = hessian(sdp.f, x);
        % linear cost
        g = jacobian(simplify(sdp.f - x'*(H/2)*x), x);
        % linear constraint
        A = jacobian(sdp_g, x);
        % constant constraint
        b = simplify(A*x - sdp_g);
        
        % get sparsity
        conic.h = sparsity(H);
        conic.a = sparsity(A);

        % create low-level conic solver
        obj.solver = casos.package.solvers.conicInternal([name '_conic'],solver,conic,opts);

        % conic problem data as function of p
        data = casadi.Function('P',{x p},{H g A b});

        % SDP problem as function of p and x
        prob = casadi.Function('S',{x p},{sdp.f sdp_g});

        % build SDP problem
        buildproblem(obj,prob,data,opts);

        % construct CasADi callback
        construct(obj,name);
    end

    function s = stats(obj)
        % Return stats.
        s = obj.solver.stats;
    end

    function print_options(obj)
        % Print list of options.
        print_options(obj.solver);
    end

    function print_option(obj,name)
        % Print information about an option.
        print_option(obj.solver,name);
    end

    function has_option(obj,name)
        % Check if option "name" exists.
        has_option(obj.solver,name);
    end
end

methods (Access=protected)
    function S = copyElement(obj)
        % Use copy constructor.
        S = casos.package.solvers.SdpsolInternal(obj);
    end
end

methods (Access={?casos.package.functions.FunctionCommon, ?casos.package.functions.FunctionWrapper})
    function S = substitute(obj,idx,expr_in,expr_out)
        % Substitute a variable for expr_in -> expr_out.
        if ischar(idx)
            % map variable name to index
            idx = index_in(obj.fhan,idx);
        end

        % get current input symbols
        var_in  = sx_in(obj.fhan);
        var_out = var_in;
        % replace variable
        var_in{idx+1}  = expr_in; % CasADi uses 0-based index
        var_out{idx+1} = expr_out;

        % substitute
        fhan = casadi.Function('f',var_in,call(obj.fhan,var_out),name_in(obj.fhan),name_out(obj.fhan));
        S = casos.package.solvers.SdpsolInternal(obj,fhan);
    end
end

methods (Access=protected)
    function [sdp, opts] = dd_relax(obj, sdp, opts)
        if isfield(opts.Kx, 'ddm') && ~isempty(opts.Kx.ddm)
            % Get DD cone sizes
            n_ddm = opts.Kx.ddm;

            % update opts
            if ~isfield(opts.Kx, 'lin')
                opts.Kx.lin = 0;
            end
            opts.add_lbc = [];
            opts.add_ubc = [];

            for i=1:length(n_ddm)

                % Extract the ddm elements from sdp.x
                ddm_x = sdp.x(size(sdp.x,1)- n_ddm(i)^2+1:end);
                ddm_x = ddm_x.reshape(n_ddm(i), n_ddm(i)); % matrix of the form (n_ddm x n_ddm)
    
                % obtain diagonal elements of ddm_x
                ddm_x_diag = diag(ddm_x);
                
                % create slack variables (only corresponding to the upper
                % triangle except the diagonal elements
                ddm_s = casadi.MX.sym('s', n_ddm(i)*(n_ddm(i)+1)/2 - n_ddm(i), 1);
                
                % create symmetric matrix with zero diagonal and upper and
                % bottom triangle equal to ddm_s            
                ind_s = triu(ones(n_ddm(i), n_ddm(i)),1); % extracts the upper triangle without the diagonal
                ind_d = ind_s + ind_s';
                
                % by the order of ddm_s, subs the 1 in ind_s with ddm_s
                s_matrix = casadi.MX(n_ddm(i)^2, 1);
    
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
    
                % new sdp.x
                sdp.x = [sdp.x; ddm_s];
    
                % new sdp.g
                n_new_g = size(sdp.g,1);
                sdp.g = [sdp.g; dd_g1(:); dd_g2(:); dd_g3(:); dd_g4(:); dd_g5(:)];
                n_new_g = size(sdp.g,1) - n_new_g;
                
                % define size of new linear cones
                opts.Kx.lin = opts.Kx.lin + n_ddm(i)^2 + size(ddm_s,1); 
                opts.Kc.lin = opts.Kc.lin + n_new_g;

                % add new bounds
                % lower bounds is always zero
                opts.add_lbc = [opts.add_lbc; zeros(n_new_g,1)];
                % upper bounds is always infinity
                opts.add_ubc = [opts.add_ubc; inf(n_new_g,1)];
            end

            % remove ddm field
            opts.Kx = rmfield(opts.Kx, 'ddm');

        end

    end

end

end
