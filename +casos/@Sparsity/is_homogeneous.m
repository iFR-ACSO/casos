function tf = is_homogeneous(obj,deg)
% Check if polynomial sparsity pattern is homogeneous (of given degree).

degsum = obj.degsum;

if nargin < 2
    % check for homogeneity of any degree
    deg = degsum(1);
end

% a polynomial is homogeneous if all terms are of equal degree
tf = all(degsum == deg);

end
