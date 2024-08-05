function f = sossol(varargin)
% Interface for convex sum-of-squares (SOS) solvers.
        
node = casos.package.solvers.sossolInternal(varargin{:});

f = casos.Function.create(node);

end
