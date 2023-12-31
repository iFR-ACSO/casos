function f = sdpsol(varargin)
% High-level interface for convex cone (SDP) solvers.
        
sol = casos.package.solvers.SdpsolInternal(varargin{:});

f = casos.Function(sol);

end
