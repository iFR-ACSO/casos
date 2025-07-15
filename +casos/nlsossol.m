function f = nlsossol(varargin)
% Interface for nonconvex sum-of-squares (SOS) solvers.
        
node = casos.package.solvers.nlsossolInternal(varargin{:});

f = casos.Function.create(node);

end
