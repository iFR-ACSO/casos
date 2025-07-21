classdef SoslmiRelaxation < casos.package.solvers.SosoptCommon
% Solve sum-of-squares problems via dual relaxation to LMI.

properties (Access=private)
    sdpsolver;

    lmiargs;
    lmi2sos;
end

properties (Constant, Access=protected)
    soslmi_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'sdpsol', 'Conic solver to be used.'
         'sdpsol_options', 'Options to be passed to the SDP solver.'
         'newton_solver', 'Solver used for the Newton simplification (defaults to the one used in sdpsol).'}
    ];
end

properties (SetAccess=private)
    class_name = 'SoslmiInterface';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SoslmiRelaxation.soslmi_options;
    end
end

methods (Access=protected)
    % internal evaluation
    out = eval_on_basis(obj,in);

    % build LMI relaxation
    buildproblem(obj,sos);
end

methods
    function obj = SoslmiRelaxation(name,sos,varargin)
        obj@casos.package.solvers.SosoptCommon(name,sos,varargin{:});

        % default options
        if ~isfield(obj.opts,'sdpsol'), obj.opts.sdpsol = 'sedumi'; end
        if ~isfield(obj.opts,'sdpsol_options'), obj.opts.sdpsol_options = struct; end
        if ~isfield(obj.opts,'newton_solver'), obj.opts.newton_solver = obj.opts.sdpsol; end

        % pass options to sdpsol
        if ~isfield(obj.opts.sdpsol_options,'error_on_fail')
            obj.opts.sdpsol_options.error_on_fail = obj.opts.error_on_fail;
        end

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

        % build LMI relaxation
        buildproblem(obj,sos);
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.sdpsolver.stats;
    end
end

end
