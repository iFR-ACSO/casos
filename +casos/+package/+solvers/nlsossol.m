function f = nlsossol(varargin)
% Interface for nonlinear sum-of-squares (SOS) solvers.
        
sol = casos.package.solvers.NlsossolInternal(varargin{:});

f = casos.Function(sol);

end