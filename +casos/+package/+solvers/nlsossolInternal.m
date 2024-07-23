function node = nlsossolInternal(name,solver,varargin)
% Internal interface for nonlinear sum-of-squares problems.

switch (solver)
    case 'sequential'
        node = casos.package.solvers.Sequential(name,varargin{:});

    case 'FeasRes'
        node = casos.package.solvers.FeasRes(name,varargin{:});

    otherwise
        error('No such sequential solver "%s".',solver)
end

end
