function f = nlsossol(varargin)
% Interface for sequential sum-of-squares (SOS) solvers.
        
node = casos.package.solvers.nlsossolInternal(varargin{:});

f = casos.Function.create(node);

end
