function node = sossolInternal(name,solver,sos,varargin)
% Internal interface for sum-of-squares solvers.

switch (solver)
    case 'alfonso'
        node = casos.package.solvers.AlfonsoSosInterface(name,sos,varargin{:});

    otherwise
        node = casos.package.solvers.SossdpRelaxation(name,solver,sos,varargin{:});
end

end
