function [Anf,bnf,cnf,Knf] = splitfreevars(A,b,c,K,tol);
% function [Anf,bnf,cnf,Knf] = splitfreevars(A,b,c,K,tol);
%
% DESCRIPTION
%   The primal form of SOS optimization specifies the decision variables
%   as free variables in the primal formulation. This results in an
%   equality constraint involving free and cone constrained variables:
%         Af*xf+Ak*xk = b  with xk \in K
%   Sedumi allows free variables in its primal formulations but many cone
%   / SDP solvers do not.  A standard trick is to replace the free
%   free variables by their positive and negative parts:
%         xf = xp-xn with xp>=0 and xn>=0
%   This function reformulates the optimization problem using this
%   variable replacement.  The result is a primal-form optimization
%   with all cone-constrained variables (no free variables).
%
%   One numerical issue with this reformulation is that the same constant
%   can be added to xp and xn without affecting x.  As a result xp and
%   xn may diverge to large values.  This can be alleviated by augmenting
%   the cost function with a small penalty on each entry of xp and xn.
%   Specifically, c'*[xf;xk] is replaced by:
%        c'*[xf;xk]+tol*||c||_inf*(xp+xn)
%   For feasibility problems, c is a vector of zeros.  In this case,
%   the cost function is modified to be tol*(xp+xn).
%
% INPUTS
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See
%       Sedumi documentation for more details.
%   tol: Penalty on each positive/negative variable added to the cost
%       function to prevent divergence (Default: tol = 1e-6*||c||_infty)
%
% OUTPUTS
%   Anf,bnf,cnf,Knf: Cone optimization problem in Sedumi format. This
%       formulation is the problem specified by the inputs but with
%       free variables split into positive and negative parts.
%
% SYNTAX
%   [Anf,bnf,cnf,Knf] = splitfreevars(A,b,c,K,tol);
%
% EXAMPLE
%

% 2/28/08  PJS     Initial Coding

if nargin< 5
    tol = [];
end
if isempty(tol)
    tol = 1e-6;
end

Anf = A;
bnf = b;
cnf = c;
Knf = K;

if ~isempty(Knf.f) && Knf.f>0
    % Replace xf with [xn xp]
    Anf = [-Anf(:,1:Knf.f) Anf];
    cnf = [-cnf(1:Knf.f);  cnf];

    % Add small penalty to prevent
    % XXX: Choosing the tol is a bit tricky. It would be cleaner to add an 
    % SOCP constraint the bounds the norm of the split free variables.
    if tol>0
        cnfinf = norm(cnf,'inf');
        if cnfinf==0
            cnf(1:2*Knf.f) = cnf(1:2*Knf.f)+tol;
        else
            cnf(1:2*Knf.f) = cnf(1:2*Knf.f)+tol*cnfinf;
        end
    end
    Knf.l = Knf.l + 2 * Knf.f;
    Knf.f = 0;
end

