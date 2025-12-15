function f = conic(varargin)
% Low-level interface for conic (SDP) solvers.

try
    sol = casos.package.solvers.conicInternal(varargin{:});
    
    f = casos.Function(sol);

catch e
    throwAsCaller(e)
end

end
