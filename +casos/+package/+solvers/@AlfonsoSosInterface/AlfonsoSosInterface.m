classdef AlfonsoSosInterface < casos.package.solvers.SosoptCommon
% Solve sum-of-squares problems using Alfonso.

properties (Access=private)
    fhan;
    ghan;
    cone;

    info = struct;
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end

properties (Constant, Access=protected)
    alfonso_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'alfonso', 'Options to be passed to Alfonso.'
         'newton_simplify', 'Perform monomial basis simplification.'}
    ];

    allow_eval_on_basis = true;
end

properties (SetAccess=private)
    class_name = 'AlfonsoInterface';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.AlfonsoSosInterface.alfonso_options;
    end
end

methods (Access=protected)
    % internal evaluation
    out = eval_on_basis(obj,in);

    % build Alfonso problem
    buildproblem(obj,sos);
end

methods
    function obj = AlfonsoSosInterface(name,sos,varargin)
        obj@casos.package.solvers.SosoptCommon(name,sos,varargin{:});

        % default options
        if ~isfield(obj.opts,'alfonso'), obj.opts.alfonso = struct; end
        if ~isfield(obj.opts,'newton_simplify'), obj.opts.newton_simplify = true; end

                % states
        if ~isfield(sos,'x')
            error('No decision variables given.')
        else
            sos.x = casos.PS(sos.x);
        end
        % parameter
        if ~isfield(sos,'p')
            % not parametrized
            sos.p = casos.PS;
        else
            sos.p = casos.PS(sos.p);
        end
        % objective
        if ~isfield(sos,'f')
            % feasibility problem
            sos.f = casos.PS(0);
        else
            sos.f = casos.PS(sos.f);
        end
        % constraints
        if ~isfield(sos,'g')
            % unconstrained optimization
            sos.g = casos.PS;
        else
            sos.g = casos.PS(sos.g);
        end

        % build Alfonso problem
        buildproblem(obj,sos);
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end
end

end
