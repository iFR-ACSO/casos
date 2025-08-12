function b = int(a,x,varargin)
% Evaluate polynomial integral operator (definite/indefinite integral).

assert(~is_operator(a), 'Not allowed for operators.')
assert(is_indet(x), 'Second argument must be vector of indeterminates.')

% number of indeterminates
nx = length(x);

switch (nargin)
    case 2
        % indefinite integral
        range = [];

    case 3
        % definite integral, range given
        range = varargin{1};

        % check size and shape of range
        % accepted: [a b], [a(:) b(:)], or [a1 ... an b1 ... bn]
        if numel(range) == 2*nx
            range = reshape(range(:)',nx,2);

        elseif length(range) == 2
            range = repmat(range(:)',nx,1);

        else
            assert(isequal(size(range),[nx 2]), 'Invalid size of range argument (got [%dx%d]).', size(range,1), size(range,2))
        end

    case 4
        % definite integral, bounds a and b given
        range = [varargin{1}(:) varargin{2}(:)];

        % check size and shape of range
        if size(range,1) == 1
            range = repmat(range,nx,1);

        else
            assert(size(range,1) == nx, 'Invalid length of upper and lower bound (got %d).', size(range,1))
        end

    otherwise
        error('Invalid number of inputs (got %d).', nargin)
end

b = a.new_poly;

% compute coefficient matrix of integral
[S,b.coeffs] = coeff_int(a.get_sparsity,a.coeffs,x,range);

% set sparsity pattern
b = set_sparsity(b,S);

end
