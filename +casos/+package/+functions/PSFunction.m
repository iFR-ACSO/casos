classdef PSFunction < casos.package.functions.FunctionInternal
% Polynomial (PS) function interface.

properties (Access=private)
    func;

    sparsity_i;
    sparsity_o;
end

properties (SetAccess=private)
    class_name = 'PSFunction';
end

methods
    function obj = PSFunction(arg1, ex_i, ex_o, name_i, name_o, varargin)
        % Create new casadi function object.
        if isa(arg1,'casadi.Function')
            % internal constructor
            name = arg1.name;

            % set function
            func = arg1;
            % set sparsity patterns
            sparsity_i = ex_i(:);
            sparsity_o = ex_o(:);

        else
            name = arg1;
        
            % parse polynomial expressions
            [sym_i,sparsity_i] = cellfun(@parse_expr, ex_i(:), 'UniformOutput', false);
            [sym_o,sparsity_o] = cellfun(@parse_expr, ex_o(:), 'UniformOutput', false);

            % define function between coefficients
            func = casadi.Function(name, sym_i, sym_o, name_i, name_o, varargin{:});
        end

        obj@casos.package.functions.FunctionInternal(name);

        obj.sparsity_i = sparsity_i;
        obj.sparsity_o = sparsity_o;
        obj.func = func;
    end

    %% Implement FunctionInternal
    function n = get_n_in(obj)
        % Number of inputs.
        n = n_in(obj.func);
    end

    function str = get_name_in(obj,varargin)
        % Name of inputs.
        str = name_in(obj.func,varargin{:});
    end

    function z = get_sparsity_in(obj,idx)
        % Monomials of inputs.
        z = obj.sparsity_i{idx+1};
    end

    function idx = get_index_in(obj,str)
        % Index of inputs.
        idx = index_in(obj.func,str);
    end

    function n = get_n_out(obj)
        % Number of outputs.
        n = n_out(obj.func);
    end

    function str = get_name_out(obj,varargin)
        % Name of outputs.
        str = name_out(obj.func,varargin{:});
    end

    function z = get_sparsity_out(obj,idx)
        % Monomials of outputs.
        z = obj.sparsity_o{idx+1};
    end

    function idx = get_index_out(obj,str)
        % Index of outputs.
        idx = index_out(obj.func,str);
    end

    function J = jacobian(obj)
        % Return Jacobian function.
        J_func = jacobian(obj.func);

        % collect input sparsity patterns for Jacobian function
        % inputs are equal to function's inputs and outputs
        sp_in = [obj.sparsity_i; obj.sparsity_o];

        % collect output sparsity pattens for Jacobian function
        % outputs are equal to function's outputs differentiated w.r.t.
        % function's inputs, in the order 
        %   jac_o0_i0,...,jac_o0_iN, ..., jac_oM_i0,...,jac_oM_iN
        sp_out = cell(length(obj.sparsity_i),length(obj.sparsity_o));

        for i=1:size(sp_out,1)
            for j=1:size(sp_out,2)
                % linear index into Jacobian's outputs
                k = sub2ind(size(sp_out),i,j);

                % build operator sparsity pattern
                sp_out{i,j} = casos.Sparsity(...
                    sparsity_out(J_func,k-1), ...
                    obj.sparsity_i{i}, ...
                    obj.sparsity_o{j} ...
                );
            end
        end

        % return Function object
        J = casos.Function.create(casos.package.functions.PSFunction(J_func,sp_in,sp_out));
    end
end

methods (Access=protected)
    %% Internal evaluation
    function argout = eval_on_basis(obj, argin)
        % Evaluate function on nonzero coordinates.
        argout = call(obj.func, argin);
    end
end

end

function [coeffs,z] = parse_expr(p)
% Return coefficients, monomials, and size of polynomial expression.

    p = casos.PS(p);
    [coeffs,z] = coordinates(p);
end
