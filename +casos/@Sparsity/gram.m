function [Z,K,Mp,Md] = gram(w)
% Return a Gram basis for a given monomial sparsity pattern.

[degmat,Lw] = get_degmat(w);

[Z,K,Mp,Md] = gram_internal(Lw,degmat,w.indets);

end
