classdef SdpsolInternal < casos.package.solvers.SolverCallback & matlab.mixin.Copyable
% Internal interface for convex cone (SDP) solvers.

properties (Constant,Access=protected)
    matrix_cones = casos.package.Cones([
        casos.package.Cones.DD;
        casos.package.Cones.SDD
    ]);
end

properties (Access=private)
    solver;
    map;    % maps from the relaxed sdp and the orignal sdp
end

properties (Access=protected)
    fhan;
    ghan;
end

methods (Access=private)
    buildproblem(obj,prob,data,opts,args);
end

methods (Static)
    function cones = get_cones
        % Return supported cones.
        cones = [
            casos.package.solvers.ConicSolver.get_cones
            casos.package.solvers.SdpsolInternal.matrix_cones
        ];
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
        end

        % ensure cone has default value
        if ~isfield(opts,'Kx')
            opts.Kx.lin = numel(sdp.x);
        elseif ~isfield(opts.Kx,'lin')
            opts.Kx.lin = 0; 
        end
        if ~isfield(opts,'Kc')
            opts.Kc.lin = numel(sdp.g);
        elseif ~isfield(opts.Kc,'lin')
            opts.Kc.lin = 0;
        end

        % relax problem to smaller easier cones (LP and SOCP)
        args = struct;
        args.dd_ubg = []; 
        args.dd_lbg = []; 
        args.dd_ubx = []; 
        args.dd_lbx = []; 
        obj.map = [];

        % rebuild problem from SDD to SOCP
        if isfield(opts.Kx,'sdd') || isfield(opts.Kc,'sdd')
            [sdp, args, map, opts] = sdd_reduce(obj, sdp, opts, args);
            map_sdd_to_primal = map;
            map_to_dual    = map.lam_x;
        else
            % if no SDD was present create an identity map
            map_sdd_to_primal.x = speye(length(sdp.x));
            map_sdd_to_primal.g = speye(length(sdp.g));
            map_to_dual     = [speye(length(sdp.x)), zeros(length(sdp.x), length(sdp.g))];
        end

        % rebuild problem from DD to LP
        if isfield(opts.Kx,'dd') || isfield(opts.Kc,'dd')
            [sdp, args, map, opts] = dd_reduce(obj, sdp, opts, args);
            map_dd_to_primal = map;
            map_to_dual   = map.lam_x;
        else
            % if no DD was present create an identity map
            map_dd_to_primal.x = speye(length(sdp.x));
            map_dd_to_primal.g = speye(length(sdp.g));
            map_to_dual    = [speye(length(sdp.x)), zeros(length(sdp.x), length(sdp.g))];
        end

        % Create the full map
        obj.map.x = map_sdd_to_primal.x*map_dd_to_primal.x;
        obj.map.g = map_sdd_to_primal.g*map_dd_to_primal.g;
        obj.map.lam_x = map_to_dual;
        obj.map.lam_a = [zeros(size(obj.map.g,1), size(obj.map.x,2)), speye(size(obj.map.g))];

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

        if isfield(sdp,'derivatives')
            % use pre-computed derivatives (undocumented)
            H = obj.map.x'*sdp.derivatives.Hf*obj.map.x;
            g = sdp.derivatives.Jf*obj.map.x;
            A = [sdp.derivatives.Jg*obj.map.x; 
                 jacobian(sdp_g(size(sdp.derivatives.Jg,1)+1:end), x)];
        else
            % quadratic cost
            H = hessian(sdp.f, x);
            % linear cost
            g = jacobian(sdp.f, x);
            % linear constraint
            A = jacobian(sdp_g, x);
        end
        % constant constraint
        b = -sdp_g;
        
        % get sparsity
        conic.h = sparsity(H);
        conic.a = sparsity(A);

        % create low-level conic solver
        obj.solver = casos.package.solvers.conicInternal([name '_conic'],solver,conic,opts);

        % conic problem data as function of p
        % NOTE: When evaluating, we must set x = 0.
        data = casadi.Function('P',{x p},{H g A b});

        % SDP problem as function of p and x
        prob = casadi.Function('S',{x p},{sdp.f sdp_g});

        % build SDP problem
        buildproblem(obj,prob,data,opts,args);

        % construct CasADi callback
        construct(obj,name);
    end

    function s = stats(obj)
        % Return stats.
        s = obj.solver.stats;
    end

    %% Options & Cones
    function print_options(obj)
        % Print list of options.
        print_options(obj.solver);
        % also print matrix cones
        print_matrix_cones(obj);
    end

    function print_option(obj,name)
        % Print information about an option.
        names = split(name,'.');

        if length(names) > 1 && ismember(names{1},{'Kx' 'Kc'}) && has(obj.matrix_cones,names{2})
            % print option
            print_one(obj.solver.get_options,names{1});
            % print matrix cone
            print_one(obj.matrix_cones,names{2});
        else
            % print option & cones
            print_option(obj.solver,name);

            if isscalar(names) && ismember(names{1},{'Kx' 'Kc'})
                % print matrix cones
                print_matrix_cones(obj);
            end
        end
    end

    function tf = has_option(obj,name)
        % Check if option "name" exists.
        tf = has_option(obj.solver,name);
    end

    function print_cones(obj)
        % Print list of supported cones.
        print_cones(obj.solver);
        % also print matrix cones
        print_matrix_cones(obj);
    end

    function print_cone(obj,name)
        % Print information about a cone.
        if has(obj.matrix_cones,name)
            % print matrix cone
            print_one(obj.matrix_cones,name);
        else
            % print cone from solver
            print_cone(obj.solver,name);
        end
    end

    function tf = has_cone(obj,name)
        % Check if cone "name" is supported.
        tf = has(obj.matrix_cones,name) || has_cone(obj.solver,name);
    end
end

methods (Access=protected)
    % reduce DD constraints to LPs
    [sdp,args,map,opts] = dd_reduce(obj,sdp,opts,args);

    % reduce SDD constraints to SOCP
    [sdp,args,map,opts] = sdd_reduce(obj,sdp,opts,args);

    % replace DD cones (constraint or decision variable form)
    [sdp,args,M_out,num_nlin,dd_index,opts] = replaceDDcones(obj,sdp,sizes,Mlin,args,opts,field)

    % replace DD cones (constraint or decision variable form)
    [sdp,args,M_out,num_nlin,sdd_index,opts] = replaceSDDcones(obj,sdp,sizes,Mlin,args,opts,field)

    function S = copyElement(obj)
        % Use copy constructor.
        S = casos.package.solvers.SdpsolInternal(obj);
    end

    function print_matrix_cones(obj)
        % Print list of supported matrix cones.
        disp('Supported Matrix Cones:')
        print_all(obj.matrix_cones);
    end
end

end
