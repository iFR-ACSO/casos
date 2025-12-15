function f = sdpsol(varargin)
% High-level interface for convex cone (SDP) solvers.
        
try
    sol = casos.package.solvers.SdpsolInternal(varargin{:});
    
    f = casos.Function(sol);

catch e
    throwAsCaller(e)
end

end
