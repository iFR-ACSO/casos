classdef SolverCallback < casadi.Callback
% Common Callback interface for custom solvers.

properties (Abstract, Access=protected)
    fhan;
    ghan;
end

methods
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

methods (Access=protected, Static)
    function [options,opts] = parse_options(options,opts)
        % Parse options.
        if nargin < 2
            opts = struct;
        end
        for fn=intersect(fieldnames(options),fieldnames(opts))'
            if isempty(fn), continue; end
            options.(fn{:}) = opts.(fn{:});
            opts = rmfield(opts,fn);
        end
    end
end

end
