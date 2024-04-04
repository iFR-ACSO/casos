function [indets,I] = combine(varargin)
% Combine indeterminate variables.

indets = casos.Indeterminates;

switch (nargin)
    case 0
        % nothing to do
        I = [];
        return

    case 1
        indets.variables = varargin{1}.variables;
        I = 1:length(indets);
        return

    case 2
        allvars = [varargin{1}.variables varargin{2}.variables];

    otherwise
        allvars = cellfun(@(a) a.variables, varargin, 'UniformOutput', false);
        allvars = [allvars{:}];
end

% remove duplicates and sort
[vars,~,I] = unique(allvars);

% return
indets.variables = vars;

end
