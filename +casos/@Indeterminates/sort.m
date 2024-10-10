function [out,ia,ic] = sort(obj)
% Sort indeterminate variables alphabetically.

[vars,ia,ic] = unique(obj.variables);   % variables are already unique

% return
out = casos.Indeterminates;
out.variables = vars;
out.transp = obj.transp;

end
