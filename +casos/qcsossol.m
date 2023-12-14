function f = qcsossol(varargin)
% Interface for quasiconvex sum-of-squares (SOS) solvers.
        
sol = casos.package.solvers.QcsossolInternal(varargin{:});

f = casos.Function(sol);

end
