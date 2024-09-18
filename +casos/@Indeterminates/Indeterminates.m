classdef Indeterminates < casos.package.core.AlgebraicObject & casos.package.core.Printable
% Indeterminate variables.

properties (GetAccess=protected, SetAccess=private)
    % cell array of strings {x1,...,xN}
    variables = {};

    % transpose flag
    transp = false;
end

properties (Dependent)
    nvars;
end

methods
    %% Public constructor
    function obj = Indeterminates(varargin)
        % Create indeterminate variables.
        if nargin == 0
            % nothing to do
            return

        elseif isa(varargin{1},'casos.Indeterminates') && nargin < 2
            % copy constructor
            obj = varargin{1};
            return
        end

        % else
        var = varargin{1};
        arg = varargin(2:end);

        % pvar / mpvar syntax
        if nargin == 1 && ischar(var)
            % syntax Indeterminates('x')
            obj.variables = {var};

        elseif ischar([arg{:}])
            % syntax Indeterminates('x','y',...)
            obj.variables = unique(varargin,'stable');

            if length(obj.variables) < nargin
                warning('Duplicate variables removed.')
            end

        elseif length(arg) <= 2
            % syntax Indeterminates('x',m,n)
            N = numel(zeros(arg{:},1));     % number of variables
            l = floor(log10(N))+1;          % number of places
            obj.variables = compose(['%s_%0' num2str(l) 'd'],var,1:N);

        else
            error('Undefined syntax.')
        end
    end

    %% Getter
    function n = get.nvars(obj)
        % Return number of variables.
        n = length(obj.variables);
    end

    function varargout = size(obj,varargin)
        % Return size of indeterminate variables.
        [varargout{1:nargout}] = size(sparse(1,obj.nvars),varargin{:});
    end

    function tf = is_equal(obj,var)
        % Check if indeterminate variables are equal.
        tf = isequal(obj.variables,var.variables);
    end

    function tf = isrow(obj)
        % Check if indeterminate variables have been transposed.
        tf = isscalar(obj) || obj.transp;
    end

    function tf = iscolumn(obj)
        % Check if indeterminate variables have not been transposed.
        tf = isscalar(obj) || ~obj.transp;
    end

    function out = str(obj)
        % Return string representation.
        out = obj.variables;
    end
end

methods (Static)
    %% Static constructor
    function v = empty()
        % Create empty indeterminate variables.
        v = casos.Indeterminates;
    end
end

methods
    % public RedefinesParen interface
    v = cat(dim,varargin);

    % AlgebraicObject interface
    function tf = is_indet(~), tf = true; end

    function z = mpower(obj,deg)
        % Return monomial(s).
        z = power(casos.package.polynomial(obj),deg);
    end

    function z = power(obj,deg)
        % Return monomial(s).
        z = power(casos.package.polynomial(obj),deg);
    end

    function obj = transpose(obj)
        % Toggle transpose flag.
        obj.transp = ~obj.transp;
    end

    % Display
    function disp(obj)
        % Display indeterminates as matrix.
        disp_matrix(obj,'()');
    end

    % Cell-like interface
    function [tf,loc] = ismember(obj,x)
        % Check if the indeterminate variables include x.
        [tf,loc] = ismember(obj.variables,x.variables);
    end
end

methods (Access=protected)
    % protected RedefinesParen interface
    obj = parenAssign(obj,idx,varargin);
    obj = parenDelete(obj,idx);
    varargout = parenReference(obj,index);
end

methods (Access={?casos.package.core.PolynomialInterface})
    % friend class access
    [indets,ic] = combine(varargin);
    [indets,ic] = sort(obj);
end

end
