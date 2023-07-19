classdef Function < casos.package.functions.FunctionWrapper
% Multiple-input, multiple-output function between polynomials.

methods
    function f = Function(name, ex_i, ex_o, name_i, name_o, varargin)
        % Create a new function.
        if nargin < 5
            % default input/output names
            name_i = compose('i%d',1:length(ex_i));
            name_o = compose('o%d',1:length(ex_o));
        end

        import casos.package.functions.*

        % parse inputs and outputs
        [type_i,args_i] = cellfun(@FunctionWrapper.parse_argument, ex_i, name_i, 'UniformOutput', false);
        [type_o,args_o] = cellfun(@FunctionWrapper.parse_argument, ex_o, name_o, 'UniformOutput', false);

        types = [type_i{:} type_o{:}];

        % select appropriate Function class
        if all(ismember(types,["DM" "SX" "DM"]))
            % fall back to casadi Function class
            wrap = CasadiFunction(name,args_i,args_o,varargin{:});

        elseif any(types == 'PS') && ~any(types == 'MX')
            % function between polynomials with symbolic coefficients
            wrap = PSFunction(name,args_i,args_o,varargin{:});

        else
            error('No matching function for inputs %s.', types);
        end

        f@casos.package.functions.FunctionWrapper(wrap,args_i,args_o,name_i,name_o);
    end

end

end