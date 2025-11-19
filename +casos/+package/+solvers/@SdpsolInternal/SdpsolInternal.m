classdef SdpsolInternal < casos.package.solvers.SolverCallback & matlab.mixin.Copyable
% Internal interface for convex cone (SDP) solvers.

properties (Access=private)
    solver;

    sdp_info;
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

    function s = info(obj)
        % Return info.
        s = obj.sdp_info;
        s.conic = obj.solver.info;
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

end
