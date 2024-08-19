function [lz,lw,h] = sossize(n,d)
% function [lz,lw,h] = sossize(n,d);
%
% DESCRIPTION 
%   This function computes the size of an SOS LMI feasibility problem:
%      p is SOS iff there is a Q>=0 such that p = z^T*Q*z
%   z is an lz-by-1 vector of monomials and Q is an lz-by-lz symmetric
%   matrix.   This can solved as an LMI in two ways:
%     Primal Form:  Find Q>=0 such A*vec(Q) = c where A is an lw-by-lz^2
%       matrix and c is an lw-by-1 vector.  A and c are obtained by
%       equating coefficients of p=z^T*Q*z
%     Dual Form:  Find lambda in R^h such that Q0+ sum_i lambda_i*N*>=0.
%       Q0 is a particular solution to p = z^T*Q*z and the N_i are a 
%       basis of symmetric homogeneous solutions satisfying z^T*Ni*z = 0.
%   This function computes the dimensions lz, lw, and h as a function
%   of the polynomial degree and number of variables. The dimensions
%   are computed assuming p is an unstructured degree d polynomial in n
%   variables, i.e. p contains all possible monomials. The SDP dimensions 
%   may be smaller for a specific polynomial.
%
% INPUTS
%   n: Number of variables in the polynomial
%   d: Degree of polynomial. 
%
% OUTPUTS
%   lz: Length of z and dimension of Q in the Gram matrix representation.
%       If d is odd then z includes all monomials of degree <= ceil(d/2).
%   lw: Number of equality constraints in the primal LMI form. 
%   h:  Number of linearly independent, symmetric, homogenous solutions 
%      satisfying z^T N z = 0.  
%
% SYNTAX 
%   [lz,lw,h] = sossize(n,d)

% 9/26/09   PJS     Initial Coding
% 10/22/10  PJS     Fix lw for odd degrees

% Number of monomials in z
dhalf = ceil(d/2);
lz = nchoosek(n+dhalf,dhalf);

% Number of monomials in w
% z'*Q*z can be a polynomial of degree 2*dhalf so this value should be 
% used to compute lw.  If d is even this has no effect. If d is odd this 
% uses the next largest even integer.

% lw = nchoosek(n+d,d);
lw = nchoosek(n+2*dhalf,2*dhalf);

% Qdim is the number of unique entries in the symmetric matrix Q
Qdim = lz*(lz+1)/2;  
h = Qdim - lw;

