function S = cat(dim,varargin)
% Concatenate polynomial sparsity patterns along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

switch (nargin-1)
    case 0
        % nothing to concatenate
        S = casos.Sparsity;

    case 1
        % single input
        S = varargin{1};

    case 2
        % two patterns given
        S1 = casos.Sparsity(varargin{1});
        S2 = casos.Sparsity(varargin{2});

        % concatenate coefficient matrices
        S = coeff_cat(S1,S2,S1.coeffs,S2.coeffs,dim);

    otherwise
        % more than two inputs
        N = (nargin-1)/2;
        % recursion
        S1 = cat(dim,varargin{1:floor(N)});
        S2 = cat(dim,varargin{ceil(N):end});
        % concatenate
        S = cat(S1,S2);
end

end
