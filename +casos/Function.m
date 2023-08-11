classdef Function < casos.package.functions.FunctionWrapper
% Multiple-input, multiple-output function between polynomials.

methods
    function f = Function(name, ex_i, ex_o, name_i, name_o, varargin)
        % Create a new function.
        import casos.package.functions.*

        if isa(name,'casos.package.functions.FunctionWrapper')
            % copy constructor
            wrap = name.wrap;

        else
        
        if nargin < 5
            % default input/output names
            name_i = compose('i%d',1:length(ex_i));
            name_o = compose('o%d',1:length(ex_o));
        end

        % parse inputs and outputs
        [type_i,args_i] = cellfun(@FunctionWrapper.parse_argument, name_i, ex_i, 'UniformOutput', false);
        [type_o,args_o] = cellfun(@FunctionWrapper.parse_argument, name_o, ex_o, 'UniformOutput', false);

        types = [type_i{:} type_o{:}];

        % select appropriate Function class
        if all(ismember(types,["DM" "SX" "DM"]))
            % fall back to casadi Function class
            wrap = CasadiFunction(name,args_i,args_o,varargin{:});

        elseif any(types == 'PS') && ~any(types == 'MX')
            % function between polynomials with symbolic coefficients
            wrap = PSFunction(name,args_i,args_o,varargin{:});

        else
            str = compose('%s',types); str(2,1:end-1) = {', '};
            error('No matching function for inputs (%s).', [str{:}]);
        end

        end

        f@casos.package.functions.FunctionWrapper(wrap);
    end
end

end