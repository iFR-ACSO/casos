function sos = sosineq(sos,symexpr,info,info2)
% SOSINEQ --- Add a new inequality constraint f(x) >= 0
%    to an SOS program 
%
% SOSP = sosineq(SOSP,EXPR)
%
% SOSP is the sum of squares program.
% EXPR is the expression on the left hand side of the constraint, i.e., f(x).
%
% EXPR can be a column vector. In this case, several inequality
% constraints will be added simultaneously to the sum of
% squares program.
%
% SOSP = sosineq(SOSP,EXPR,[a b]) is used to specify the interval 
% a <= x <= b where f(x) has to be non-negative. Currently, this is 
% possible only if the sum of squares program is univariate, and the input
% is of type 'sym'. The lower limit a can be -Inf, and similarly b can be +Inf.
%
% SOSP = sosineq(SOSP,EXPR,'sparse') will use the fact that a polynomial 
% is sparse to compute a sum of squares decomposition using an optimally 
% reduced set of monomials, which makes the size of the resulting semidefinite 
% program smaller.
%
% SOSP = sosineq(SOSP,EXPR,'sparsemultipartite',VARCELLS) will enable
% efficient exploitation of sparse multipartite structure. In this case,
% VARCELLS is a cell which describes the partition of the independent
% variables. For example
%
%    SOSP = sosineq(SOSP,EXPR,'sparsemultipartite',{[x,y],[z]})
%
% corresponds to multipartite structure where the partition of the
% independent variables are {x,y} and {z}.
% 

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
%                                      M. Peet (6), D. Jagt (6)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
% (6) Cybernetic Systems and Controls Laboratory, Arizona State University,
%     Tempe, AZ 85287-6106, USA.
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes

% 03/20/02 - SP -- Interval info
% 12/15/21 - DJ -- Interval adjustment for dpvar case



if isfield(sos,'symvartable')
    
	% Interval
	if nargin > 2 & isnumeric(info)
        if length(sos.symvartable) ~= 1
            error('Not a univariate polynomial.');
        end;
        for i = 1:size(symexpr,1)
            symexpr(i) = preprocess(symexpr(i),sos.symvartable,info);
            sos = sosconstr(sos,'ineq',symexpr(i));
            sos.expr.type{sos.expr.num} = 'posint';
        end;
        return;
	end;

elseif nargin > 2 && isnumeric(info) % DJ, 12/15/21
    if length(sos.vartable) ~= 1
        error('Not a univariate polynomial.');
    end
    for i = 1:size(symexpr,1)
        xvar = pvar(sos.vartable{1});
        symexpr(i) = polypreprocess(symexpr(i),xvar,info);
        sos = sosconstr(sos,'ineq',symexpr(i));
        sos.expr.type{sos.expr.num} = 'posint';
    end
    return;

end;
    
% Sparse polynomial
if nargin > 2 & strcmp(info,'sparse');
   for i = 1:size(symexpr,1)
       sos = sosconstr(sos,'ineq',symexpr(i));
       sos.expr.type{sos.expr.num} = 'sparse';
   end;
   return;
end;

% Sparse polynomial
if nargin > 2 & strcmp(info,'sparsemultipartite');
   for i = 1:size(symexpr,1)
       sos = sosconstr(sos,'ineq',symexpr(i));
       sos.expr.type{sos.expr.num} = 'sparsemultipartite';
       sos.expr.multipart{sos.expr.num} = info2;
   end;
   return;
end;

sos = sosconstr(sos,'ineq',symexpr);

% ============================================================
function newsymexpr = preprocess(symexpr,var,info)
 
% Get the maximum degree of the independent variable
maxdeg = 0;
dummy = diff(symexpr,var);
while dummy ~= 0
    maxdeg = maxdeg+1;
    dummy = diff(dummy,var);
end;

% Substitute var
if info(2)==Inf
    newvar = var+info(1);
    newsymexpr = subs(symexpr,var,newvar);
elseif info(1)==-Inf
    newvar = -var+info(2);
    newsymexpr = subs(symexpr,var,newvar);
else
    newvar = (info(2)-info(1))/2*(1-var)/(1+var) + (info(2)+info(1))/2;
    newsymexpr = subs(symexpr,var,newvar)*(1+var).^maxdeg;
end;

% ============================================================
function pnewsymexpr = polypreprocess(psymexpr,var,info)
 
% Get the maximum degree of the independent variable
maxdeg = 0;
dummy = diff(psymexpr,var);
while any(dummy.C ~= 0)
    maxdeg = maxdeg+1;
    dummy = diff(dummy,var);
end

% Substitute var
if info(1)==-Inf && info(2)==Inf
    pnewsymexpr = psymexpr;
elseif info(2)==Inf
    newvar = var+info(1);
    pnewsymexpr = subs(psymexpr,var,newvar);
elseif info(1)==-Inf
    newvar = -var+info(2);
    pnewsymexpr = subs(psymexpr,var,newvar);
else
    error('Imposing inequality on finite interval is not supported for pvar/dpvar implementation')
end
