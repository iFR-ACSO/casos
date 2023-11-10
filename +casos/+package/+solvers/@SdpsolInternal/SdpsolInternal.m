classdef SdpsolInternal < casos.package.solvers.SolverCallback
% Internal interface for convex cone (SDP) solvers.

properties (Access=private)
    solver;
end

properties (Access=protected)
    fhan;
    ghan;
end

methods (Access=private)
    buildproblem(obj,prob,data);
end

methods
    argout = eval(obj,argin);

    function obj = SdpsolInternal(name,solver,sdp,opts)
        obj@casos.package.solvers.SolverCallback;

        % options
        if nargin < 4
            opts = struct;
        elseif isfield(opts,'Kg')
            % map constraint cones
            opts.Ka = opts.Kg;
            opts = rmfield(opts,'Kg');
        end

        % decision variables
        x = sdp.x;
        % parameter
        if isfield(sdp,'p')
            p = sdp.p;
        else
            p = [];
        end

        % quadratic cost
        H = hessian(sdp.f, x);
        % linear cost
        g = jacobian(simplify(sdp.f - x'*(H/2)*x), x);
        % linear constraint
        A = jacobian(sdp.g, x);
        % constant constraint
        b = simplify(A*x - sdp.g);
        
        % get sparsity
        conic.h = sparsity(H);
        conic.a = sparsity(A);

        % create low-level conic solver
        obj.solver = casos.package.solvers.conicInternal([name '_conic'],solver,conic,opts);

        % conic problem data as function of p
        data = casadi.Function('P',{p},{H g A b});

        % SDP problem as function of p and x
        prob = casadi.Function('S',{x p},{sdp.f sdp.g});

        % build SDP problem
        buildproblem(obj,prob,data);

        % construct CasADi callback
        construct(obj,name);
    end
end

end
