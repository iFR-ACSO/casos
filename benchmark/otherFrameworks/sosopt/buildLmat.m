function [L,Rbasis] = buildLmat(z);  
% function [L,Rbasis] = buildLmat(z);
%
% DESCRIPTION 
% Builds a matrix representation of the linear operator, L:
%      L: Q --> z(x)'*Q*z(x)
% The domain is the set of nxn matrices (not necessarily symmetric) and
% the range space is the set of polynomials.  In constructing L, the 
% basis used for the range space is given by the monomials in the vector 
% Rbasis.  The basis used for the domain is Q_i = reshape(e_i,n,n) where
% e_i is the i^th standard basis vector in R^{n*n}.  
% Thus the construction of L yields:
%            z(x)'*Q*z(x) =   Rbasis'*L*Q(:)
%   
% INPUTS 
%   z:  lz-by-1 vector of monomials 
%
% OUTPUTS 
%   L: lr-by-lz^2 matrix representation of the linear operator L: Q -->
%       z(x)'*Q*z(x). 
%   Rbasis: lr-by-1 vector of monomials which form a basis for the 
%       range space of L.  
%
% SYNTAX 
%   [L,Rbasis] = buildLmat(z);
%
% EXAMPLE
%   pvar x1 x2 
%   z = [x1; x2; x1^2];
%   [L,Rbasis]=buildLmat(z); 
%   Q=rand(3); 
%   z'*Q*z,  Rbasis'*L*Q(:), z'*Q*z-Rbasis'*L*Q(:)

% 1/28/08  PJS     Re-coding based on LOCALbuildL

% Error checking
z = polynomial(z);
if ~ismonom(z) || size(z,2)~=1
    error('z must be a nz-by-1 vector of monomials');
end
   
if z.nvar ==0
    % Handle case where z(x) = constant
    Rbasis = polynomial(1);
    L = 1;
    return;
else    
    % Build vector of all monomials in range of L, [z'*Q_i*z]
    % XXX This step can cause memory issues (zQz is large matrix)
    zmat = z.coef'*z.degmat;
    lz = size(zmat,1);
    zQz = repmat(zmat,lz,1) + kron(zmat,ones(lz,1));

    % The basis for the range space is the set of unique monomials in zQz
    [Rbasis,temp,idx] = unique(zQz,'rows');  % Rbasis(idx,:) = zQz
    lr = length(Rbasis);
    chkval = 0; % skip validity check
    Rbasis = polynomial(speye(lr) , Rbasis, z.var, [lr,1],chkval);

    % Build matrix representation of the linear operator, L:
    %      L: Q --> z'*Q*z
    if 1
        % Every column of L has the non-zero entry 1 in exactly one row
        Qdim = lz^2;
        L = sparse(idx,1:Qdim,1);
    else
        % PJS 12/6/09
        % Q is assumed to be symmetric.  In this case the columns of L
        % corresponding to the strict upper triangle of Q are zero, the
        % cols corresponding to the strict lower triangle have exactly
        % one entry  = 2, and the cols corresponding to the main diag of
        % Q have exactly one entry = 1.
        % 
        % This version of L is more sparse (fewer non-zero entries)
        % and I thought it would result in a speed-up by Sedumi.  
        % The results were identical on sosoptdemo3 using the 'image'
        % form (both comp time and step-by-step Sedumi print-outs of each 
        % iteration). I'll stick with the original version above because 
        % this version causes issues with the 'kernel' form.
        Qdim = lz^2;
        Lvals = sparse(lz^2,1);
        [ridx,cidx]=ind2sub([lz lz],1:Qdim);
        Lvals( ridx==cidx ) = 1;
        Lvals( ridx>cidx ) = 2;
        L = sparse(idx,1:Qdim,Lvals);                       
    end    
end
