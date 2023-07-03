function [indets,dga,dgb] = combineVar(indetA,indetB,dga,dgb)
% Combine variables of two polynomials and extend degree matrices.

nta = size(dga,1);
ntb = size(dgb,1);

% combine variables
[indets,~,ic] = unique([indetA indetB]);

nv = length(indets);

% get nonzero degrees
[ia,ja,da] = find(dga);
[ib,jb,db] = find(dgb);

% extend degree matrices to combined variables
dga = sparse(ia,ic(ja),da,nta,nv);
dgb = sparse(ib,ic(length(indetA)+jb),db,ntb,nv);

end
