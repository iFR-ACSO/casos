function f = sossol(varargin)
% Interface for convex sum-of-squares (SOS) solvers.

try
    node = casos.package.solvers.sossolInternal(varargin{:});
    
    f = casos.Function.create(node);

catch e
    throwAsCaller(e)
end

end
