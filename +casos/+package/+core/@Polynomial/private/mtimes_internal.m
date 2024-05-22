function [S,coeffs] = mtimes_internal(S1,S2,coeff1,coeff2)
% Compute coefficient matrix for matrix multiplication.

nta = S1.nterm;
ntb = S2.nterm;

sza = size(S1);
szb = size(S2);

if isrow(S1) && isvector(S2)
    % inner vector product
    % (sum_a c_a'*x^a)*(sum_b c_b*x^b) = sum_a sum_b (c_a'*c_b)*(x^a*x^b)
    coeffs = reshape(coeff1*T(coeff2),nta*ntb,1);

else
    % Vectorized code to compute coef matrix
    % from multipoly
    idx1 = reshape(1:sza(1)*sza(2),[sza(1),sza(2)]);
    idx1 = repmat(idx1,[nta 1]);
    idx1 = idx1(:);
    
    idx2 = repmat(1:nta,[sza(1) 1]);
    idx2 = repmat(idx2(:),[sza(2) 1]);
    
    idx = sub2ind([sza(1)*sza(2) nta],idx1,idx2);
    acoefcol = T(coeff1);
    acoefcol = reshape(acoefcol(idx),[nta*sza(1) sza(2)]);
    
    bcoefcol = reshape(T(coeff2),[szb(1) szb(2)*ntb]);
    
    tempcoef = acoefcol*bcoefcol;
    
    idx1 = reshape(1:sza(1)*nta,[sza(1),nta]);
    idx1 = repmat(idx1,[ntb*szb(2) 1]);
    idx1 = idx1(:);
    
    idx2 = repmat(1:(ntb*szb(2)),[sza(1) 1]);
    idx2 = repmat(idx2(:),[nta 1]);
    
    idx = sub2ind([nta*sza(1) ntb*szb(2)],idx1,idx2);
    coeffs = T(reshape(tempcoef(idx),sza(1)*szb(2),nta*ntb));
end

% update coefficients and degrees
[S,coeffs] = coeff_mtimes(S1,S2,coeffs);

end
