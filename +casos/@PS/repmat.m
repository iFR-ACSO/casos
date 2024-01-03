function obj = repmat(obj,varargin)
% Repeat copies of polynomial.
assert(nargin > 0, 'Not enough input arguments.')

% repition scheme
rep = horzcat(varargin{:});

assert(isrow(rep) && length(rep) == 2, 'Replication factors must be a pair (row) of integers or two integer scalars.')

% new dimensions
new_sz = rep.*obj.matdim;

% repeat coefficients
obj.coeffs = reshape(repmat(obj.coeffs,rep),obj.nterm,prod(new_sz));
% set dimensions
obj.matdim = new_sz;

end
