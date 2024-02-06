function f = sossol(varargin)
% Interface for convex sum-of-squares (SOS) solvers.
        
sol = casos.package.solvers.SossolInternal(varargin{:});

f = casos.Function(sol);

end
