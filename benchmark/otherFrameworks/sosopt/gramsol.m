function [z,Q0,N,C] = gramsol(p,g,z)  
% function [z,Q0,N,C] = gramsol(p,g,z);  
%
% DESCRIPTION 
%   Construct particular and null solutions for Gram matrix constraint.
%   Given a 1-by-1 polynomial p and a 1-by-N polynomial g then the 
%   following Gram matrix decomposition represents a set of linear 
%   constraints on the unknown coefficients c and the unknown entries of 
%   the matrix Q:
%           p(x)+g(x)*c = z(x)'*Q*z(x)
%   This function computes Q0, Ni and Ci such that:
%         z'*Q0*z = p, z'*Ni*z=0, and z'*Ci*z = gi
%   Thus all solutions of the Gram matrix decomposition are given by:
%         Q0 + sum_i lambda_i N_i + sum_j c_j C_j
%
% INPUTS 
%   p: 1-by-1 polynomial
%   g: 1-by-N polynomial [Optional]
%   z: lz-by-1 monomials vector used in Gram matrix constraint [Optional]
%
% OUTPUTS 
%   z: lz-by-1 monomials vector used in Gram matrix constraint
%   Q0: lz-by-lz matrix satisfying z'*Q0*z = p
%   N: lz^2-by-M matrix whose i^th column satisfies 
%          z'*reshape(N(:,i),[lz lz])*z = 0.
%   C: lz^2-by-N matrix whose i^th column satisfies
%          z'*reshape(C(:,i),[lz lz])*z = gi.
%
% SYNTAX 
%   [z,Q0,N] = gramsol(p);  
%   [z,Q0,N,C] = gramsol(p,g);  
%   [z,Q0,N,C] = gramsol(p,g,z);  
%
% EXAMPLE
%   pvar x1 x2 
%   p = 5*x1+3*x1*x2-4*x2^2*x1;
%   g = [x1^2-9, 6*x2*x1];
%   [z,Q0,N,C] = gramsol(p,g);  
%
%   lz = length(z);
%   p - z'*Q0*z
%   for i1=1:size(N,2)
%       z'*reshape(N(:,i1),lz,lz)*z
%   end
%   for i1=1:size(C,2)
%       g(i1) - z'*reshape(C(:,i1),lz,lz)*z
%   end

% 2/29/08  PJS Initial Coding
% 11/7/10  PJS Rewrote to call GRAMCONSTRAINT to construct matrix eqns   


%----------------------------------------------------------------------
% p(x)+g(x)*c = z(x)'*Q*z(x) is equivalent to:
%   Rbasis'*(p_R + g_R*c) = Rbasis'*L*Q(:)
%----------------------------------------------------------------------
if nargin < 2
    g = [];
    [z,p_R,L,Rbasis] = gramconstraint(p);
elseif nargin<3
    [z,p_R,L,g_R,Rbasis] = gramconstraint(p,g);
else
    [z,p_R,L,g_R,Rbasis] = gramconstraint(p,g,z);
end
lz = length(z);
lr = length(Rbasis);

%----------------------------------------------------------------------
% Find symmetric Q0 and Ci such that 
%   p(x) = z(x)'*Q0*z(x) and gi(x) = z(x)'*Ci*z(x)
%----------------------------------------------------------------------
Q0 = L \ p_R;
Q0 = LOCALSymmetrize(Q0);
Q0 = reshape(Q0,[lz lz]);
if ~isempty(g)
    C = L \ g_R;
    C = LOCALSymmetrize(C);
else
    C = [];
end

%----------------------------------------------------------------------
% Build a basis for nullspace of L.
% The columns of N contain vec(Ni) where z'*Ni*z = 0	  
%----------------------------------------------------------------------
N = sparse( lz^2 , 0);
for i1=1:lr
    idx = find(L(i1,:));
    lidx = length(idx);
    if lidx>1
        % Null space vectors can be found row-by-row due to structure of L
        % XXX The non-zero entries of L are all 1 which means L(i1,idx)
        % must always be of the form ones(1,N).  N1compact can then
        % be explicitly given by [ -ones(1,N-1); eye(N-1)]
        %N1compact = null(full(L(i1,idx)),'r');
        N1compact = [-ones(1,lidx-1); eye(lidx-1)];
        N1 = sparse( lz^2 , size(N1compact,2) );
        N1(idx,:) = N1compact;
        
        % Symmetrize Null Space Solutions.  Anti-symmetric solutions are
        % symmetrized to all-zeros and removed. Then remove non-unique 
        % solutions.
        % XXX I think this can be done faster and w/o the call to unique
        N1 = LOCALSymmetrize(N1);
        idx = find( sum(abs(N1))==0 );  
        N1(:,idx) = [];
        N1 = unique(N1','rows')';        
        N = [N N1];
    end
end

%----------------------------------------------------------------------
% function Mout = LOCALSymmetrize(Min);
%   Min is an n^2-by-m matrix with the k^th column representing a
%   vectorized n-by-n matrix, Mk.  Mout is an n^2-by-m matrix with 
%   the k^th column representing the symmetric form of Mk.
%----------------------------------------------------------------------
function Mout = LOCALSymmetrize(Min)

[nsq,m]=size(Min);
n = sqrt(nsq);
Mout = sparse(nsq,m);
for i1=1:m
    Mk = reshape(Min(:,i1),n,n);
    Mk = (Mk+Mk')/2;
    Mout(:,i1) = Mk(:);
end
