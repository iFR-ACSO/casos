function p = cat(dim,varargin)
% Concatenate polynomials along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

switch (nargin-1)
    case 0
        % nothing to concatenate
        p = [];

    case 1
        % single input
        p = varargin{1};

    case 2
        % two patterns given
        p1 = varargin{1};
        p2 = varargin{2};

        p = p1.new_poly;
        
        % concatenate coefficient matrices
        [S,p.coeffs] = coeff_cat(get_sparsity(p1),get_sparsity(p2),p1.coeffs,p2.coeffs,dim);
        
        % set sparsity
        p = set_sparsity(p,S);

    otherwise
        % more than two inputs
        N = (nargin-1)/2;
        % recursion
        p1 = cat(dim,varargin{1:floor(N)});
        p2 = cat(dim,varargin{ceil(N):end}); 
        % concatenate
        p = cat(dim,p1,p2);
end

end
