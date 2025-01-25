function sol = conicInternal(name,solver,conic,varargin)
% Internal interface for conic (SDP) solvers.
        
switch lower(solver)
    case 'sedumi'
        % solve conic problem using SeDuMi
        sol = casos.package.solvers.SedumiInterface(name,conic,varargin{:});

    case 'mosek'
        % solve conic problem using MOSEK
        sol = casos.package.solvers.MosekInterface(name,conic,varargin{:});

    case 'scs'
        % solve conic problem using SCS
        sol = casos.package.solvers.SCSInterface(name,conic,varargin{:});

    case 'clarabel'
        % solve conic problem using SCS
        sol = casos.package.solvers.ClarabelInterface(name,conic,varargin{:});
        
    otherwise
        error('Solver "%s" undefined.', solver)
end

end
