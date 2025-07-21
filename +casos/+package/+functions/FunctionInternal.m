classdef (Abstract) FunctionInternal < casos.package.functions.FunctionCommon
% Internal interface for callable functions.

properties (Abstract,SetAccess=private)
    class_name;
end

properties (SetAccess = private)
    name;
end

methods (Abstract)
    % inputs
    n = get_n_in(obj);
    s = get_name_in(obj,i);
    z = get_sparsity_in(obj,i);
    i = get_index_in(obj,str);

    % outputs
    n = get_n_out(obj);
    s = get_name_out(obj,i);
    z = get_sparsity_out(obj,i);
    i = get_index_out(obj,str);
end

methods
    function obj = FunctionInternal(name,varargin)
        % Superclass constructor.
        obj@casos.package.functions.FunctionCommon(varargin{:});
        obj.name = name;
    end

    %% Getter
    function z = get_monomials_in(obj,i)
        % Return monomials of input arguments.
        warning('DEPRECATED. Use sparsity instead.')
        z = monomials(get_sparsity_in(obj,i));
    end

    function val = get_default_in(obj,i) %#ok<INUSD>
        % Return default values of input arguments.
        val = 0;
    end

    function sz = get_size_in(obj,i)
        % Return size of input arguments.
        sz = size(get_sparsity_in(obj,i));
    end

    function z = get_monomials_out(obj,i)
        % Return monomials of output arguments.
        warning('DEPRECATED. Use sparsity instead.')
        z = monomials(get_sparsity_out(obj,i));
    end

    function sz = get_size_out(obj,i)
        % Return size of output arguments.
        sz = size(get_sparsity_out(obj,i));
    end

    function s = get_stats(~)
        % Return empty stats.
        s = struct;
    end

    function J = jacobian(obj) %#ok<STOUT>
        % Return Jacobian function if supported.
        error('Derivatives cannot be calculated for %s.',obj.name)
    end

    %% Call internal
    function argout = call(obj,argin,on_basis)
        % Call function.
        if nargin > 2 && on_basis
            % evaluate on coordinates
            argout = eval_on_basis(obj,argin);

        else
            % get nonzero coordinates for basis
            in = cellfun(@(p,i) poly2basis(p,get_sparsity_in(obj,i)), argin(:), num2cell((1:get_n_in(obj))'-1), 'UniformOutput', false);
            % evaluate on coordinates
            out = eval_on_basis(obj,in);
            % return polynomials as result
            argout = cellfun(@(c,i) casos.package.polynomial(get_sparsity_out(obj,i),c), out(:), num2cell((1:get_n_out(obj))'-1), 'UniformOutput', false);
        end
    end
end

methods (Access=protected)
    %% Internal evaluation
    function argout = eval_on_basis(obj,argin) %#ok<STOUT,INUSD>
        % Evaluate function on nonzero polynomials.
        error('Not implemented.')
    end
end

end
