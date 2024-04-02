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
        if length(range) == 2
            range = repmat(range(:)',nx,1);

        elseif numel(range) == 2*nx
            range = reshape(range(:)',nx,2);

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

% combine variables
[indets,dgf,~] = combineVar(f.indets,x.indets,f.degmat,x.degmat);

% find position if variables
[tf,xloc] = ismember(indets,x.indets);

% integral(sum_a c_a*x1^a1...xn^an) 
%   = sum_a c_a*integral(x^1...xn^an)
%   = sum_a c_a/((a1+1)...(an+1)) x1^(a1+1)...xn^(an+1)
g = casos.PS;

% increase degrees
dgf(:,tf) = dgf(:,tf) + 1;

% compute scaling term for each monomial
scale = prod(dgf(:,tf),2);

% scale coefficients
coeffs = scale.\f.coeffs;

if isempty(range)
    % indefinite integral
    degmat = dgf;

else
    % definite integral
    % integral(c_A*x^A, r1, r2) = c_A*(r2^A - r1^A)
    deg0 = dgf(:,tf)';

    % select range expressions
    r = x.coeffs(xloc(tf),:)*range;

    % compute exponents
    R1 = sx_prod(r(:,1).^deg0,1);
    R2 = sx_prod(r(:,2).^deg0,1);

    % coefficients of sum_A c_A*(r2^A - r1^A)
    coeffs = coeffs.*(R2 - R1)';
    % degree matrix
    degmat = dgf(:,~tf);

    % select remaining variables
    indets = indets(~tf);

    % make degree matrix unique
    [coeffs,degmat] = uniqueDeg(coeffs, degmat);

    % remove zero terms
    [coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);
end

% new polynomial
g.coeffs = coeffs;
g.degmat = degmat;
g.indets = indets;
g.matdim = f.matdim;

end
