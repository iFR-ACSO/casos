classdef SossdpRelaxation < casos.package.solvers.SosoptCommon
% Solve sum-of-squares problems by relaxation to SDP.

properties (Access=private)
    sdpsolver;

    gram_x;
    gram_g;

    basis_x_out;
    basis_g_out;
end

properties (Constant,Access=protected)
    sossdp_options = [casos.package.solvers.SosoptCommon.sosopt_options
        {'sdpsol_options', 'Options to be passed to the SDP solver.'}
    ];
end

properties (SetAccess=protected)
    class_name = 'SossdpRelaxation';
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SossdpRelaxation.sossdp_options;
    end
end

methods
    out = call(obj,in);

    function obj = SossdpRelaxation(name,solver,sos,varargin)
        obj@casos.package.solvers.SosoptCommon(name,sos,varargin{:});

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

        % build SDP problem
        buildproblem(obj,solver,sos);
    end

    function s = get_stats(obj)
        % Return stats.
        s = obj.sdpsolver.stats;
    end
end

end
