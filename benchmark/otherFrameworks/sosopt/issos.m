function [feas,z,Q,f] = issos(p,opts)
% function [feas,z,Q,f] = issos(p,opts)
%
% DESCRIPTION 
%   This function tests if a polynomial p is a sum of squares. If p is SOS
%   then the function returns the Gram matrix decomposition, z^TQz with
%   Q>=0. The function also returns an SOS decomposition, sum_i f(i)^2.
%
% INPUTS 
%   p:  1-by-1 polynomial
%   opts: Options for optimization.  See SOSOPTIONS for more details.
%       
% OUTPUTS
%   feas: feas=1 if the polynomial is a SOS and feas=0 otherwise.
%   z: lz-by-1 vector of monomials in the Gram matrix decomposition of p, 
%        if it exists.  If feas=0 then z =[].
%   Q: lz-by-lz positive semidefinite matrix in the Gram matrix 
%        decomposition of p, if it exists. If feas=0 then Q=[].
%   f: lf-by-1 vector of polynomials in the SOS decomposition of p, if
%        it exists.  In this case p = f(1)^2 + f(2)^2 + ... + f(end)^2.
%        If feas=0 then f=[].
%
% SYNTAX 
%   [feas,z,Q,f] = issos(p)
%   [feas,z,Q,f] = issos(p,opts)
%
% EXAMPLE 
%  pvar x1 x2;
%  p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;
%  [feas,z,Q,f] = issos(p)
%  
%  % Verify Gram matrix decomposition
%  p - z'*Q*z
%  % Verify Q >=0
%  min(eig(Q))
%  % Verify SOS decomposition of p
%  p - f'*f
%
% See also: sosopt, sosoptions

% 6/3/09   PJS     Initial Coding

% Setup
if nargin==1
    opts = sosoptions;
end
szp = size(p);
if ~all( szp==[1 1] )
    error('Input polynomial p must be 1-by-1.');
end

% Convert polyconstr to a one-side polynomial
if isa(p,'polyconstr')
    if strcmp(p.RelOp,'==')
        error('polconstr inputs with RelOp set to == are not allowed.');
    else
        p = p.OneSide;
    end
end

% Promote to poly and handle constant case
p = polynomial(p);
if p.nvar==0
   p = double(p);
   if p>=0
        feas = 1;
        z = polynomial(1);
        Q = p;
        f = polynomial(sqrt(p));
        return       
   else
        feas = 0;
        z = [];
        Q = [];
        f = [];
        return
   end
end

% Call sosopt and unpack
[info,dopt,sossol] = sosopt( p>=0 , opts);
feas = info.feas;
if feas==1
    z = sossol(1).z;
    Q = sossol(1).Q;
    
    % Construct SOS decomposition 
    % Q = U*T*U' --> Q = L'*L where L = sqrt(T)*U'
    [U,T] = schur(full(Q));  
    T = max(T,0);  % This should only be zeroing tiny negative evals
    L = sqrt(T)*U';
    f = L*z; 
else
    z = [];
    Q = [];
    f = [];
end
