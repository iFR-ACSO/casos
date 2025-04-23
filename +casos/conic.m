function f = conic(varargin)
% Low-level interface for conic (SDP) solvers.
        
sol = casos.package.solvers.conicInternal(varargin{:});

f = sol;
%f = casadi.Function(sol);

end
