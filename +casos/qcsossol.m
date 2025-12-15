function f = qcsossol(varargin)
% Interface for quasiconvex sum-of-squares (SOS) solvers.
        
try
    node = casos.package.solvers.qcsossolInternal(varargin{:});
    
    f = casos.Function.create(node);

catch e
    throwAsCaller(e)
end

end
