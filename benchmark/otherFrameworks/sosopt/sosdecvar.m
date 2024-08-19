function [p,D,z] = sosdecvar(dstr,z)
% function [s,D,z] = sosdecvar(dstr,z)
%
% DESCRIPTION 
%   This function creates a SOS decision variable, s, for use in SOS
%   optimization problems.  s has the form s = z'*D*z where D is a 
%   symmetric matrix of decision variable coefficients. 
%   
%   Note: s is not automatically enforced to be SOS by SOSOPT.  The list
%   of constraints passed to SOSOPT must include s>=0.
%
% INPUTS 
%   dstr: Character string used in creating the decision variables.
%   z: Nz-by-1 column vector of monomials.
%
% OUTPUTS
%   s: 1-by-1 polynomial.
%   D: Nz-by-Nz symmetric matrix of pvars. The pvar names are specified by
%        dstr, e.g. if dstr='d' then the (i,j) entry of D is d_i_j for j>i.
%   z: Nz-by-1 column vector of monomials. This is the same as the input z.
%
% SYNTAX
%   s = sosdecvar(dstr,z)
%   [s,D,z] = sosdecvar(dstr,z)
%
% See also: polydecvar, sosopt, gsosopt


%  PJS 4/25/2009  Initial Coding
%  PJS 5/29/2009  Modification for a polydecvar class
%  PJS 10/28/2010 Elimination of polydecvar class. Created SOSDECVAR fcn.


if nargin == 0
    p = polynomial;
    D = polynomial(zeros(0,1));
    z = polynomial(zeros(0,1));    
    return
elseif nargin==2
    type = 'vec';
elseif nargin~=3
    errstr1 = 'Invalid syntax for the "polydecvar" command.';
    errstr2 = ' Type "help polydecvar" for more information.';
    error([errstr1 errstr2]);
end

z = polynomial(z);
if ~isa(z,'polynomial') || ~ismonom(z)   
    error('Second input must be a vector of monomials')
end

% Create SOSDECVAR using Gram Matrix form

% Get z degmat with same ordering as the monomials appear in z 
z = polynomial(z(:));
zdeg = z.coefficient'*z.degmat;
lz = length(z);

% Create (symmetric) coefficient matrix
D = mpvar(dstr,lz,lz,'s');
%Dvar = char(D(:));
Dvar = char(D);
Dvar = Dvar(:);

% Create polynomial and combine terms due to symmetry
pcoef = sparse( ones(lz^2,1) );

nxvar = size(zdeg,2);
pdeg = repmat(zdeg,[lz 1]) + ...
    reshape(repmat(zdeg(:)',[lz 1]),[lz^2 nxvar]);    
pdeg = [speye(lz^2) pdeg];

pvarname = [Dvar(:); z.var];
p = polynomial(pcoef,pdeg,pvarname,[1 1]);
p = combine(p);
