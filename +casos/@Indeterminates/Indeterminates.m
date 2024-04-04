classdef Indeterminates < casos.package.core.AlgebraicObject
% Indeterminate variables.

properties (GetAccess=protected, SetAccess=private)
    % cell array of strings {x1,...,xN}
    variables = {};
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
            N = numel(zeros(arg{:}));
            obj.variables = compose('%s_%d',var,1:N);

        else
            error('Undefined syntax.')
        end
    end

    %% Getter
    function varargout = size(obj,varargin)
        % Return size of indeterminate variables.
        [varargout{1:nargout}] = size(sparse(1,length(obj.variables)),varargin{:});
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
end

methods (Access=protected)
    % protected RedefinesParen interface
    obj = parenAssign(obj,idx,varargin);
    obj = parenDelete(obj,idx);
    varargout = parenReference(obj,index);
end

methods (Access={?casos.package.core.AlgebraicObject})
    % friend class access
    [indets,ic] = combine(varargin);
    [indets,ic] = sort(obj);
end

end
