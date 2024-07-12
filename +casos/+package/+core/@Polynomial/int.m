function g = int(f,x,varargin)
% Evaluate polynomial integral operator (definite/indefinite integral).

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

g = f.new_poly;

% compute coefficient matrix of integral
[S,g.coeffs] = coeff_int(f.get_sparsity,f.coeffs,x,range);

% set sparsity pattern
g = set_sparsity(g,S);

end
