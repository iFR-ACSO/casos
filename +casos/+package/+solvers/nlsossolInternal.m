function node = nlsossolInternal(name,solver,varargin)
% Internal interface for noncovex sum-of-squares problems.

switch (solver)
    case 'filter'
        node = casos.package.solvers.FilterLinesearch(name,varargin{:});
    otherwise
        error('No such nonlinear solver "%s".',solver)
end

end
