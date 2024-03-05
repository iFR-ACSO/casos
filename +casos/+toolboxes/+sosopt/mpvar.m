function p = mpvar(str,varargin)
% Create a vector or matrix of indeterminate variables.

if nargin > 1 && ischar(varargin{end})
    warning('Options are not supported.')
    varargin(end) = [];
end

% new polynomial
p = casos.PS(str,varargin{:});

end
