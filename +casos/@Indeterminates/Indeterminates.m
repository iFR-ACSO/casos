% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (InferiorClasses = {?casadi.DM, ?casadi.SX, ?casadi.MX}) ...
    Indeterminates < casos.package.core.AlgebraicObject & casos.package.core.Printable
% Indeterminate variables.
%
%% Constructor summary:
%
%   Indeterminates()
%
% create empty list of indeterminate variables.
%
%   Indeterminates(char, ...)
%   Indeterminates(char, int)
%
% create list of indeterminate variables.
%
%% Static constructor summary:
%
%   [x1, ...] = create(char, ...)
%   [x1, ...] = create(char, int)
%
% create indeterminate variables individually.
%
%%

properties (GetAccess=protected, SetAccess=private)
    % cell array of strings {x1,...,xN}
    variables = cell(1,0);

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
            variables = {var};

        elseif ischar([arg{:}])
            % syntax Indeterminates('x','y',...)
            variables = unique(varargin,'stable');

            if length(variables) < nargin
                warning('Duplicate variables removed.')
            end

        elseif length(arg) <= 2
            % syntax Indeterminates('x',m,n)
            N = numel(zeros(arg{:},1));     % number of variables
            l = floor(log10(N))+1;          % number of places
            variables = compose(['%s_%0' num2str(l) 'd'],var,1:N);

        else
            error('Undefined syntax.')
        end

        % ensure that names are stored as row vector
        obj.variables = reshape(variables,1,[]);
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

    function varargout = create(varargin)
        % Create indeterminate variables individually.
        %
        %   [x1, ...] = create(char, ...)
        %   [x1, ...] = create(char, int)
        %
        varargout = tuple2cell(casos.Indeterminates(varargin{:}));
    end
end

methods
    % public RedefinesParen interface
    v = cat(dim,varargin);

    % AlgebraicObject interface
    function tf = is_indet(~), tf = true; end

    function z = mpower(obj,deg)
        % Return sum of monomial(s).
        %
        % For an even degree p, this corresponds to the p-norm ||x||_p 
        % of the vector x to the power of p.
        if ~isscalar(deg)
            error('Power of indeterminates must be scalar.'); 
        end
        
        % else:
        z = sum(power(obj,deg));
    end

    function z = power(obj,deg)
        % Return monomial(s).
        z = power(casos.package.polynomial(obj),deg);
    end

    function z = prod(obj,~)
        % Return product of variables.
        z = prod(casos.package.polynomial(obj));
    end

    function c = subs(a,x,b)
        % Substitute indeterminate variables.
        c = subs(casos.package.polynomial(a),x,casos.package.polynomial(b));
    end

    function z = sum(obj,~)
        % Return sum of variables.
        z = sum(casos.package.polynomial(obj));
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

    % Convert to cell of single variables
    function cell_of_indets = tuple2cell(obj)
        % Return a cell array of individual indeterminate variables.
        cell_of_indets = cellfun(@(var) casos.Indeterminates(var), obj.variables, 'UniformOutput', false);
    end
end

methods (Access={?casos.package.core.PolynomialInterface})
    % friend class access
    [indets,ic] = combine(varargin);
    [indets,ia,ic] = sort(obj);
end

end
