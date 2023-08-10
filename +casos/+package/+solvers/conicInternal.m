function sol = conicInternal(name,solver,conic,varargin)
% Internal interface for conic (SDP) solvers.
        
switch (solver)
    case 'sedumi'
        % solve conic problem using SeDuMi
        sol = casos.package.solvers.SedumiInterface(name,conic,varargin{:});

    otherwise
        error('Solver "%s" undefined.', solver)
end

end
