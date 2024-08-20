function [z,b,Aq,Ac,Rbasis] = gramconstraint(p,g,z);  
% function [z,b,Aq,Ac,R] = gramconstraint(p,g,z);  
%
% DESCRIPTION 
%   Create matrices which can be used to enforce the Gram matrix
%   decomposition of a polynomial.  Given a 1-by-1 polynomial p and
%   a 1-by-N polynomial g then the following Gram matrix decomposition
%   represents a set of linear constraints on the unknown coefficients
%   c and the unknown entries of the matrix Q:
%           p(x)+g(x)*c = z(x)'*Q*z(x)
%   This function computes monomial vectors z and R, matrics Ac and Aq,
%   and a vector b such that
%         p(x) = R(x)'*b  
%         g(x) = R(x)'*Ac
%         z(x)'*Q*z(x) = R(x)'*Aq*Q(:)
%   Thus the polynomial constraint p(x)+g(x)*c = z(x)'*Q*z(x) can
%   be equivalently expressed as the matrix equation:
%         b + Ac*c = Aq*Q(:)  
%
% INPUTS 
%   p: 1-by-1 polynomial
%   g: 1-by-N polynomial [Optional]
%   z: lz-by-1 monomials vector used in Gram matrix constraint [Optional]
%
% OUTPUTS 
%   z: monomial vector used in Gram matrix constraint
%   b, Aq, Ac: b is a vector and Aq, Ac are two matrices such that the
%     matrix equation Aq*Q(:) + Ac*c = b is equivalent to the polynomial
%     constraint equation p(x)+g(x)*c = z(x)'*Q*z(x)
%   R: lr-by-1 basis of monomials. 
%
% SYNTAX 
%   [z,b,Aq,R] = gramconstraint(p);  
%   [z,b,Aq,Ac,R] = gramconstraint(p,g);  
%
% EXAMPLE
%   pvar x1 x2 
%   p = 5*x1+3*x1*x2-4*x2^2*x1;
%   g = [x1^2-9, 6*x2*x1];
%   [z,b,Aq,Ac,R] = gramconstraint(p,g);  
%
%   lz = length(z);
%   p-R'*b
%   g-R'*Ac
%   Q = [0 2.5 0 1.5; 2.5 0 0 0; 0 0 0 -2; 1.5 0 -2 0];
%   z'*Q*z - R'*Aq*Q(:)

% 1/28/08  PJS     Initial Coding


if nargin < 2
    g = [];
    z = [];
elseif nargin<3
    z = [];
end

%----------------------------------------------------------------------
% Step 1) Build monomials vector needed to write constraint in terms
%   of the the Gram matrix: 
%           p(x) + g(x)*c = z(x)'*Q*z(x)
%----------------------------------------------------------------------
if ~isempty(z)
    lz = length(z);
else
    h = [p g];
    mxdg = max(h.maxdeg)/2;
    mndg = min(h.mindeg)/2;
    z = monomials( h.nvar , [floor(mndg):ceil(mxdg)]);

    % Discard monomials based on simple checks
    % (Further reduction is performed in SOSSIMPLIFY)
    lz = size(z,1);
    mxdg = ceil(max(h.degmat,[],1)/2);
    [ridx,cidx] = find( repmat(mxdg,lz,1) < z );
    z(ridx,:) = [];

    lz = size(z,1);
    mndg = floor(min(h.degmat,[],1)/2);
    [ridx,cidx] = find( z< repmat(mndg,lz,1) );
    z(ridx,:) = [];

    lz = size(z,1);
    chkval = 0; % skip validity check
    z = polynomial(speye(lz),z,h.var,[lz 1],chkval);
end
if nargout==1
    return
end

%----------------------------------------------------------------------
% Step 2) Build a matrix representation, L, of the linear operator 
%   which maps matrices to polynomials, Q --> z'*Q*z.  L and
%   Rbasis satisfy z'*Q*z=Rbasis'*L*Q(:)
%----------------------------------------------------------------------
[L,Rbasis] = buildLmat(z);  
lr = length(Rbasis);

%----------------------------------------------------------------------
% Step 3) Write p(x) and g(x) in terms of Rbasis:  
%   p(x) = Rbasis'*p_R and g(x) = Rbasis'*g_R
%----------------------------------------------------------------------
p_R = poly2basis(p,Rbasis);
if ~isempty(g)
    g_R = poly2basis(g,Rbasis);
else
    g_R = [];
end

%----------------------------------------------------------------------
% Step 4) p(x)+g(x)*c = z(x)'*Q*z(x) is equivalent to:
%   Rbasis'*(p_R + g_R*c) = Rbasis'*L*Q(:)
% This is equivalent to the matrix/vector relation
%   L*Q(:) = p_R + g_R*c 
% Rewrite this in the form:
%   Aq*Q(:) = b + Ac*c
%----------------------------------------------------------------------
b = p_R;
Aq = L;
Ac = g_R;
if nargin==1
    Ac = Rbasis;    
end
