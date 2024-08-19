function Asym = symconstraints(A,K)
% function Asym = symconstraints(A,K);
%
% DESCRIPTION
%   Sedumi primal form has equality constraints of the form
%        A*x=b
%   where x is a vector of free and cone constrained decision variables
%   described by K. Sedumi allows for the equality constraint to represent
%   non-symmetric constraints on the entries of positive semidefinite
%   matrix decision variables.  SOSOPT converts some polynomial
%   constraints into non-symmetric equality constraints. Some solvers,
%   e.g. SDPT3 and SDPLR, expect symmetric equality constraints.  This
%   function symmetrizes the equality constraints on the decision
%   variables that are restricted to the positive semidefinite cone.
%
% INPUTS
%   A,K:  Cone optimization data in Sedumi format.  See
%       Sedumi documentation for more details.
%
% OUTPUTS
%   Asym: Cone optimization data representing symmetric equality
%       constraints on the positive semidefinite decision variables.
%
% SYNTAX
%   Asym = symconstraints(A,K);
%
% EXAMPLE
%

% 11/16/10  PJS     Initial Coding

ptr = 0;
if isfield(K,'f') && ~isempty(K.f)
    ptr = ptr + K.f;
end
if isfield(K,'l') && ~isempty(K.l)
    ptr = ptr + K.l;
end
Ns = length(K.s);

Asym = A;
for i1=1:Ns
    % Find single indices into strict lower and upper half of matrix
    n = K.s(i1);
    [ridx,cidx] = find(tril(ones(n),-1)); 
    lidx = ptr+sub2ind([n n],ridx,cidx);
    uidx = ptr+sub2ind([n n],cidx,ridx);
    
    % Symmetrize Constraints
    tmp = ( Asym(:,lidx) + Asym(:,uidx) ) / 2;    
    Asym(:,lidx) = tmp;
    Asym(:,uidx) = tmp;
    
    % Update pointer
    ptr = ptr+n^2;
end

