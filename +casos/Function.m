classdef Function < casos.package.functions.FunctionWrapper
% Multiple-input, multiple-output function between polynomials.

methods
    function f = Function(name, ex_i, ex_o, name_i, name_o, varargin)
        % Create a new function.
        import casos.package.functions.*

        try
        if nargin < 1
            % null function
            wrap = [];

        elseif isa(name,'casos.package.functions.FunctionWrapper')
            % copy constructor
            wrap = name.wrap;

        elseif isa(name,'casadi.Function')
            % wrap casadi Function object
            wrap = CasadiFunction(name);

        else
        
        if nargin < 5
            % default input/output names
            name_i = compose('i%d',0:length(ex_i)-1);
            name_o = compose('o%d',0:length(ex_o)-1);
        end

        % parse inputs and outputs
        type_i = cellfun(@FunctionWrapper.parse_argument, ex_i);
        type_o = cellfun(@FunctionWrapper.parse_argument, ex_o);

        types = [type_i type_o];

        % select appropriate Function class
        if all(ismember(types,["DM" "SX" "MX"]))
            % fall back to casadi Function class
            wrap = CasadiFunction(name,ex_i,ex_o,name_i,name_o,varargin{:});

        elseif ~any(types == 'MX')
            % function between polynomials with symbolic coefficients
            wrap = PSFunction(name,ex_i,ex_o,name_i,name_o,varargin{:});

        else
            str = compose('%s',types); str(2,1:end-1) = {', '};
            error('No matching function for inputs (%s).', [str{:}]);
        end

        end

        catch e
            throwAsCaller(e)
        end

        f@casos.package.functions.FunctionWrapper(wrap);
    end
end

methods (Static)
    function f = create(node)
        % Create a new function.
        f = set_wrapped(casos.Function, node);
    end
end

end
