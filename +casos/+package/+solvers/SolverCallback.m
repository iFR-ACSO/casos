classdef SolverCallback < casadi.Callback
% Common Callback interface for custom solvers.

properties (Abstract, Access=protected)
    fhan;
    ghan;
end

properties (Constant, Access=protected)
    solver_options = casos.package.Options({'error_on_fail', 'Throw exceptions when function evaluation fails (default true).'});
end

properties (Access=protected)
    opts = struct;
end

methods (Static, Access=protected)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SolverCallback.solver_options;
    end
end

methods
    function obj = SolverCallback(opts)
        % Instantiate base class.
        if nargin > 0
            % assert that all options exist
            check(obj.get_options, opts);
            % store options
            obj.opts = opts;
        end

        % default options
        if ~isfield(obj.opts,'error_on_fail')
            obj.opts.error_on_fail = true;
        end
    end

    %% Common Callback interface
    function n = get_n_in(obj)
        % Return number of inputs.
        n = n_in(obj.fhan);
    end

    function n = get_n_out(obj)
        % Return number of outputs.
        n = n_out(obj.ghan);
    end

    function str = get_name_in(obj,i)
        % Return names of input arguments.
        str = name_in(obj.fhan,i);
    end

    function str = get_name_out(obj,i)
        % Return names of output arguments.
        str = name_out(obj.ghan,i);
    end

    function sp = get_sparsity_in(obj,i)
        % Return sparsity of input arguments.
        sp = sparsity_in(obj.fhan,i);
    end

    function sp = get_sparsity_out(obj,i)
        % Return sparsity of output arguments.
        sp = sparsity_out(obj.ghan,i);
    end
end

end
