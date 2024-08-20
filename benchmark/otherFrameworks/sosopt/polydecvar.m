function [p,D,w] = polydecvar(dstr,w,type)
% function [p,D,w] = polydecvar(dstr,w)
%
% DESCRIPTION 
%   This function creates a polynomial decision variable, p, for use in
%   SOS optimization problems. p has the form p = D'*w where D is a column 
%   vector of decision variable coefficients.
%
% INPUTS 
%   dstr: Character string used in creating the decision variables.
%   w: Nw-by-1 column vector of monomials.
%
% OUTPUTS
%   p: 1-by-1 polynomial.
%   D: D is an Nw-by-1 column vector of pvars with names specified by 
%        dstr, e.g. if dstr = 'd' then the i^th entry of D is d_i.
%   w: Nw-by-1 column vector of monomials. This is the same as the input w.
%
% SYNTAX
%   p = polydecvar(dstr,w)
%   [p,D,w] = polydecvar(dstr,w)
%
% See also: sosdecvar, sosopt, gsosopt

%  PJS 4/25/2009  Initial Coding
%  PJS 5/29/2009  Modification for a polydecvar class
%  PJS 10/28/2010 Elimination of polydecvar class. Created SOSDECVAR fcn.

if nargin == 0
    p = polynomial;
    D = polynomial(zeros(0,1));
    w = polynomial(zeros(0,1));    
    return
elseif nargin==2
    type = 'vec';
elseif nargin~=3
    errstr1 = 'Invalid syntax for the "polydecvar" command.';
    errstr2 = ' Type "help polydecvar" for more information.';
    error([errstr1 errstr2]);
end

w = polynomial(w);
if ~isa(w,'polynomial') || ~ismonom(w)   
    error('Second input must be a vector of monomials')
end

if strcmp(type,'vec');
    % Create poly dec var using coefficient vector form

    % Get w degmat with same ordering as the monomials appear in w 
    w = polynomial(w(:));
    wdeg = w.coefficient'*w.degmat;
    lw = length(w);
            
    % Create coefficient vector
    D = mpvar(dstr,lw,1);
    Dvar = D.varname;
        
    % Create polynomial
    pcoef = sparse(ones(lw,1));
    pdeg = [speye(lw) wdeg];
    pvarname = [Dvar; w.var];
    chkval = 0; % skip validity check
    p = polynomial(pcoef,pdeg,pvarname,[1 1],chkval);        
    p = combine(p);
elseif strcmp(type,'mat');
    % Create poly dec var using Gram Matrix form
    % 10/28/10: This is now implemented in SOSDECVAR function.
    % TYPE is undocumented in POLYDECVAR but retained in code
    % for backwards compatibility.
    [p,D,w] = sosdecvar(dstr,w);
else
    error('Type must be either ''vec'' or ''mat''.')
end
