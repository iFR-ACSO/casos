function f = qcsossol(varargin)
% Interface for quasiconvex sum-of-squares (SOS) solvers.
        
node = casos.package.solvers.qcsossolInternal(varargin{:});

f = casos.Function.create(node);

end
