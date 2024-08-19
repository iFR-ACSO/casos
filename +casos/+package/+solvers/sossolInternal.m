function node = sossolInternal(name,solver,sos,varargin)
% Internal interface for sum-of-squares solvers.

node = casos.package.solvers.SossdpRelaxation(name,solver,sos,varargin{:});

end
