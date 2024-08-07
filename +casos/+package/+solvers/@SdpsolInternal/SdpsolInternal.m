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

        if isfield(sdp,'derivatives')
            % use pre-computed derivatives (undocumented)
            H = sdp.derivatives.Hf;
            g = sdp.derivatives.Jf;
            A = sdp.derivatives.Jg;

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

        % NOTE: As per Casadi v3.6.5, functions' input and output names 
        % must be mutually exclusive unless explicity permitted.
        fopt = struct('allow_duplicate_io_names',true);

        % substitute
        fhan_new = casadi.Function('f',var_in,call(obj.fhan,var_out),name_in(obj.fhan),name_out(obj.fhan),fopt);
        S = casos.package.solvers.SdpsolInternal(obj,fhan_new);
    end
end

end
