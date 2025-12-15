function f = nlsossol(varargin)
% Interface for nonconvex sum-of-squares (SOS) solvers.

try
    node = casos.package.solvers.nlsossolInternal(varargin{:});
    
    f = casos.Function.create(node);

catch e
    throwAsCaller(e)
end

end
