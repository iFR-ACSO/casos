function node = qcsossolInternal(name,solver,varargin)
% Internal interface for quasiconvex sum-of-squares problems.

switch (solver)
    case 'bisection'
        node = casos.package.solvers.QuasiconvBisection(name,varargin{:});

    otherwise
        error('No such quasiconvex solver "%s".',solver)
end

end
