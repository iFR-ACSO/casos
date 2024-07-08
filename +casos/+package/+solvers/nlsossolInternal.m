function node = nlsossolInternal(name,solver,varargin)
% Internal interface for quasiconvex sum-of-squares problems.

switch (solver)
    case 'sequential'
        node = casos.package.solvers.Sequential(name,varargin{:});

    otherwise
        error('No such sequential solver "%s".',solver)
end

end
