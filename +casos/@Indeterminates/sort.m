function [out,I] = sort(obj)
% Sort indeterminate variables alphabetically.

[vars,I] = sort(obj.variables);

% return
out = casos.Indeterminates;
out.variables = vars;

end
