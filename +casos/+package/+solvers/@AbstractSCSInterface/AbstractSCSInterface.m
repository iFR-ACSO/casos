classdef (Abstract) AbstractSCSInterface < casos.package.solvers.ConicSolver
% Abstract interface for SCS-like conic solvers.

properties (Access=protected)
    fhan;
    ghan;
    cone;

    info;
end

methods (Abstract, Static, Access=protected)
    % return triangular sparsity pattern
    S = sparsity_triangular(k);
end

methods (Abstract, Access=protected)
    % call conic solver
    [x,y,s] = call_solver(obj,data,K);
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = AbstractSCSInterface(name,conic,varargin)
        % Construct conic interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});
    end

    function s = stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end
end

methods (Access=protected)
    buildproblem(obj);

    %% Quasi-Static helper functions
    % SCS-style (de)-vectorization of the semidefinite cone
    V = sdp_vec(obj,M,varargin);
    M = sdp_mat(obj,V,varargin);
end

end
