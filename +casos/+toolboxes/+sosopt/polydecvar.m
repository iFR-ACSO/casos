function p = polydecvar(dstr,w,varargin)
% Create a polynomial decision variable.

if nargin > 2
    warning('Type argument is not supported.')
end

% new symbolic polynomial
p = casos.PS.sym(dstr,w);

end
