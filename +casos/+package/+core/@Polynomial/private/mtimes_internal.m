function [S,coeffs] = mtimes_internal(S1,S2,coeff1,coeff2)
% Compute coefficient matrix for matrix multiplication.

nta = S1.nterm;
ntb = S2.nterm;

sza = size(S1);     % MxL
szb = size(S2);     % LxN

% (sum_a C_a*x^a)*(sum_b C_b*x^b) = sum_a sum_b (C_a*C_b)*(x^a*x^b)
%
% coefficient matrices are stored as (see casos.Sparsity)
%
%        | a11 ... aM1 ... a1L ... aML |
%   C1 = |  :       :       :       :  |
%        | c11 ... cM1 ... c1L ... cML |
%
% and
%
%        | x11 ... xL1 ... x1N ... xLN |
%   C2 = |  :       :       :       :  |
%        | z11 ... zL1 ... z1N ... zLN |
%
% The core idea is to reshape C1 and reshape the transpose of C2 
% to multiply
%
%   | a11 ... a1L |
%   |  :       :  |
%   | c11 ... c1L |   | x11 ... x1N ... z11 ... z1N |
%   |  :       :  | * |  :       :       :       :  |
%   | aM1 ... aML |   | xL1 ... xLN ... zL1 ... zLN |
%   |  :       :  |
%   | cM1 ... cML |
%
% to obtain
%
%        | a1'*x1 ... a1'*xN ... a1'*z1 ... a1'*zN |
%        |    :          :          :          :   |
%        | c1'*x1 ... c1'*xN ... c1'*z1 ... c1'*zN |
%   T1 = |    :          :          :          :   |
%        | aM'*x1 ... aM'*xN ... aM'*z1 ... aM'*zN |
%        |    :          :          :          :   |
%        | cM'*x1 ... cM'*xN ... cM'*z1 ... cM'*zN |
%
% then reshape to obtain
%
%        | a1'*x1 ... a1'*z1 |
%        |    :          :   |
%        | c1'*x1 ... c1'*z1 |
%        |    :          :   |
%        | aM'*x1 ... aM'*z1 |
%        |    :          :   |
%        | cM'*x1 ... cM'*z1 |
%   T2 = |    :          :   |
%        | a1'*xN ... a1'*z1 |
%        |    :          :   |
%        | c1'*xN ... c1'*z1 |
%        |    :          :   |
%        | aM'*xN ... aM'*z1 |
%        |    :          :   |
%        | cM'*xN ... cM'*z1 |
%
% lastly transpose and reshape to
%
%       | a1'*x1 ... aM'*x1 ... a1'*xN ... aM'*xN |
%       |    :          :          :          :   |
%   C = | a1'*z1 ... aM'*z1 ... a1'*zN ... aM'*zN |
%       |    :          :          :          :   |
%       | c1'*x1 ... cM'*x1 ... c1'*xN ... cM'*xN |
%       |    :          :          :          :   |
%       | c1'*z1 ... cM'*z1 ... c1'*zN ... cM'*zN |
%
% in short:
%
%   C = R((R(C1)*R(C2'))')
%

if isrow(S1) && iscolumn(S2)
    % inner vector product (M = N = 1)
    % no need to reshape C1 and C2, hence
    %
    %   C = R((C1*C2')') = R(C2*C1')
    %
    temp = coeff2*T(coeff1);
    coeffs = reshape(temp,nta*ntb,1);

elseif iscolumn(S2)
    % matrix-vector product (N = 1)
    % no need to reshape C2, hence
    %
    %   C = R((R(C1)*C2')') = R(C2*R(C1)')
    % 
    temp = coeff2*reshape(coeff1,sza(1)*nta,sza(2))';
    coeffs = reshape(temp,nta*ntb,sza(1));

else
    % matrix-matrix product
    tempa = reshape(coeff1,sza(1)*nta,sza(2));
    tempb = reshape(T(coeff2),szb(1),szb(2)*ntb);
    tempc = reshape(tempa*tempb,sza(1)*nta*szb(2),ntb);
    coeffs = reshape(tempc',nta*ntb,sza(1)*szb(2));

% else
%     % Vectorized code to compute coef matrix
%     % from multipoly
%     idx1 = reshape(1:sza(1)*sza(2),[sza(1),sza(2)]);
%     idx1 = repmat(idx1,[nta 1]);
%     idx1 = idx1(:);
% 
%     idx2 = repmat(1:nta,[sza(1) 1]);
%     idx2 = repmat(idx2(:),[sza(2) 1]);
% 
%     idx = sub2ind([sza(1)*sza(2) nta],idx1,idx2);
%     acoefcol = T(coeff1);
%     acoefcol = reshape(acoefcol(idx),[nta*sza(1) sza(2)]);
% 
%     bcoefcol = reshape(T(coeff2),[szb(1) szb(2)*ntb]);
% 
%     tempcoef = acoefcol*bcoefcol;
% 
%     idx1 = reshape(1:sza(1)*nta,[sza(1),nta]);
%     idx1 = repmat(idx1,[ntb*szb(2) 1]);
%     idx1 = idx1(:);
% 
%     idx2 = repmat(1:(ntb*szb(2)),[sza(1) 1]);
%     idx2 = repmat(idx2(:),[nta 1]);
% 
%     idx = sub2ind([nta*sza(1) ntb*szb(2)],idx1,idx2);
%     coeffs = T(reshape(tempcoef(idx),sza(1)*szb(2),nta*ntb));
end

% update coefficients and degrees
[S,coeffs] = coeff_mtimes(S1,S2,coeffs);

end
