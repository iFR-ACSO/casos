function [indets,dg1,dg2,dg3,dg4] = combineVars(S1,S2,S3,S4)
% Combine variables of four sparsity patterns and extend degree matrices.

nvi = cumsum([0 S1.nvars S2.nvars S3.nvars]);

% combine variables
[indets,ic] = combine(S1.indets,S2.indets,S3.indets,S4.indets);

nv = length(indets);

% get nonzero degrees
[i1,j1,d1] = find(S1.degmat);
[i2,j2,d2] = find(S2.degmat);
[i3,j3,d3] = find(S3.degmat);
[i4,j4,d4] = find(S4.degmat);

% extend degree matrices to combined variables
dg1 = sparse(i1,ic(nvi(1)+j1),d1,S1.nterm,nv);
dg2 = sparse(i2,ic(nvi(2)+j2),d2,S2.nterm,nv);
dg3 = sparse(i3,ic(nvi(3)+j3),d3,S3.nterm,nv);
dg4 = sparse(i4,ic(nvi(4)+j4),d4,S4.nterm,nv);

end
