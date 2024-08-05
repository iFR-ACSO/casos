function tf = is_wellposed(obj)
% Check if sparsity is well posed.

tf = size(obj.coeffs,2) == prod(obj.matdim)      ... number of elements
    && size(obj.coeffs,1) == size(obj.degmat,1)  ... number of terms
    && size(obj.degmat,2) == length(obj.indets); ... number of variables

end
