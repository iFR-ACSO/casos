function [GAM,vars,xopt] = findbound(p,ineq,eq,DEG,options)
% FINDBOUND --- Find a global/constrained lower bound for a polynomial. 
%
% [GAM,VARS,XOPT] = findbound(P,OPTIONS)
%
% P is a polynomial of even degree, whose global lower bound is to be
% computed. The function computes the largest GAM such that (P - GAM) is a
% sum of squares.
%
% OPTIONS is an optional argument which specifies the solver and the
% solver-specific parameters on its fields as
%   options.solver
%   options.params
% the default value for options.solver is 'SeDuMi'. The solver parameters
% are fields of options.params as, for instance, options.params.tol = 1e-9.
%
%
% In addition, a vector VARS containing all the variables in P, and the
% optimal argument XOPT (in the same ordering as VARS) will  be computed,
% such that by substituting VARS = XOPT into P  we obtain P = GAM. If no
% such optimal argument is found, the function  will return an empty XOPT.
%
% [GAM,VARS,XOPT] = findbound(P,INEQ,EQ,DEG,OPTIONS) will compute a lower bound 
% for the constrained optimization problem: minimize P subject to INEQ >= 0
% and EQ = 0. Here INEQ and EQ are cells of polynomials which define the
% inequality and equality constraints. The function computes a lower bound
% based on Schmudgen's positivstellensatz, i.e., it computes the largest
% GAM such that
% 
%   (P-GAM) - LAM'*EQ - SIGMA1'*INEQ - ....
% 
% is a sum of squares, where LAM is a vector of polynomials, and SIGMA is a
% vector of sums of squares. The degree of the expression will be
% determined by the input argument DEG. Higher value of DEG will return a
% better lower bound, although the computational cost will also increase.
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

% switched gam to dpvar -MMP 8/6/2021
% Adjusted to use symvar instead of findsym - DJ, 10/14/2021
% Removed num2cell in case when variables are pvar - DJ, 12/10/2021
% Fix for polynomial with empty degmat + dpvar implementation - DJ, 04/24/2022

switch nargin
    case 1 
        options.solver='sedumi';
        ineq = [];
        eq = [];
    case 2 
        options = ineq;
        ineq = [];
        eq = [];
    case 3
        options.solver='sedumi';
        ineq = ineq(:);
        DEG = eq;
        eq = [];
    case 4
        if ~isstruct(DEG)
            ineq = ineq(:);
            eq = eq(:);
            options.solver='sedumi';
        else
            options = DEG;
            ineq = ineq(:);
            DEG=eq;
            eq = [];
        end
    case 5 
        ineq = ineq(:);
        eq = eq(:);
end

% Inequality and equality expressions should not contain any decision
% variables
if isa(ineq,'dpvar')
    ineq = compress(ineq);
    if ~isempty(ineq.dvarname)
        warning('The inequality constraints should not contain any decision variables; converting to independent variables...')
    end
    ineq = dpvar2poly(ineq);
end
if isa(eq,'dpvar')
    eq = compress(eq);
    if ~isempty(eq.dvarname)
        warning('The equality constraints should not contain any decision variables; converting to independent variables...')
    end
    eq = dpvar2poly(eq);
end
% Decision variables in the expression should not appear in the equality
% and inequality constraints
if isa(p,'dpvar')
    if (isa(ineq,'polynomial') && any(ismember(p.dvarname,ineq.varname))) || ...
            (isa(eq,'polynomial') && any(ismember(p.dvarname,eq.varname)))
        warning('Certain independent variables in the constraints appear as decision variables in the polynomial; converting to independent variables')
        p = dpvar2poly(p);
    end
end
    
vect = [p; ineq; eq];

% Find the independent variables, check the degree
if isa(vect,'sym')
   %varschar = findsym(vect);  %vars = sym(['[',varschar,']']);
   vars = symvar(vect);         % DJ, 10-14-2021
   nvars = size(vars,2);
   if nargin > 2
       degree = 2*floor(DEG/2);
       deg = zeros(length(vect),1);
       for i = 1:length(vect)
           deg(i) = double(feval(symengine,'degree',vect(i),converttochar(vars)));
           if deg(i) > degree
               error('One of the expressions has degree greater than DEG.');
           end;
       end;   
   else
       % Can change, to call maple only once.
       for var = 1:nvars;
           if rem(double(feval(symengine,'degree',p)),2) ;
               disp(['Degree in ' char(vars(var)) ' should be even. Otherwise the polynomial is unbounded.']);
               GAM = -Inf;
               xopt = [];
               return;
           end;
       end;
   end;
   syms gam;
   dvars = gam;
   
elseif isa(vect,'polynomial')
   varname = vect.varname;
   vars = [];
   for i = 1:size(varname,1)
       pvar(varname{i});
       vars = [vars eval(varname{i})];
   end;
   
   if nargin > 2
       degree = 2*floor(DEG/2);
       degmat = sum(vect.degmat,2);
       deg = zeros(length(vect),1);
       for i = 1:length(vect)
           idx = vect.coefficient(:,i)~=0;
           deg_i = max(degmat(idx));
           if isempty(deg_i)    % DJ, 04/24/2022
               deg(i) = 0;
           else
               deg(i) = deg_i;
           end
           if deg(i) > degree
               error('One of the expressions has degree greater than DEG.');
           end;
       end;   
   else
       deg = mod(max(p.degmat,[],1),2);
       if sum(deg)~=0
           i = find(deg~=0);
           disp(['Degree in ' varname{i(1)} ' should be even. Otherwise the polynomial is unbounded.']);
           GAM = -Inf;
           xopt = [];
           return;
       end;
   end;
   dpvar gam;
   dvars = gam;
   
elseif isa(vect,'dpvar')
   varname = vect.varname;
   vars = combine(polynomial(varname'));
   dvarname = vect.dvarname;
   dvars = combine(dpvar(dvarname'),'extended');
   ndvars = length(dvarname);
   
   if nargin > 2
       degree = 2*floor(DEG/2);
       degmat = sum(vect.degmat,2);
       deg = zeros(length(vect),1);
       for i = 1:length(vect)
           i_indcs = ((i-1)*ndvars + 1 : i*ndvars);
           idx = vect.C(i_indcs,:)~=0;
           deg_i = max(degmat(idx));
           if isempty(deg_i)
               deg(i) = 0;
           else
               deg(i) = deg_i;
           end
           if deg(i) > degree
               error('One of the expressions has degree greater than DEG.');
           end
       end 
   else
       deg = mod(max(p.degmat,[],1),2);
       if sum(deg)~=0
           i = find(deg~=0);
           disp(['Degree in ' varname{i(1)} ' should be even. Otherwise the polynomial is unbounded.']);
           GAM = -Inf;
           xopt = [];
           return;
       end
   end
   dpvar gam;
   dvars = [dvars,gam];
else
    error('Function and constraints should be specified as "sym", "polynomial" or "dpvar" class objects');
end

% Construct other valid inequalities
if length(ineq)>1
    for i = 1:2^length(ineq)-1
        Ttemp = dec2bin(i,length(ineq));
        T(i,:) = str2num(Ttemp(:))';
    end;
    T = T(find(T*deg(2:1+length(ineq)) <= degree),:);
    
    deg = [deg(1); T*deg(2:1+length(ineq)); deg(2+length(ineq):end)];
    for i = 1:size(T,1)
        ineqtempvect = (ineq.').^T(i,:);
        ineqtemp(i) = ineqtempvect(1);
        for j = 2:length(ineqtempvect)
            ineqtemp(i) = ineqtemp(i)*ineqtempvect(j);
        end;
    end;
    ineq = ineqtemp;
end;

prog = sosprogram(vars,dvars);
expr = p-gam;
for i = 1:length(ineq)
    [prog,sos] = sossosvar(prog,monomials(vars,0:floor((degree-deg(i+1))/2)));
    expr = expr - sos*ineq(i);
end;
for i = 1:length(eq)
    [prog,pol] = sospolyvar(prog,monomials(vars,0:degree-deg(i+1+length(ineq))));
    expr = expr - pol*eq(i);
end;
prog = sosineq(prog,expr);
prog = sossetobj(prog,-gam);
[prog,info] = sossolve(prog,options);


xopt = [];

if (info.dinf>1e-2) || (info.pinf>1e-2)
    disp('No lower bound could be computed. Unbounded below or infeasible?');
    GAM = '-inf';
    return;
else
    GAM = double(sosgetsol(prog,gam,16));
    % Returning the solution
    % The code here will only work in the rank one case.
    % Otherwise, we'll solve the f_i(x)=0 (efficiently)
    %

    % Find where the variables are (perhaps they're in the wrong order? Ask Stephen.
    ix = find(sum(prog.extravar.Z{1},2)==1) ;
    xopt = prog.solinfo.extravar.dual{1}(ix,1) ;
    vars = vars.';
        
    if isa(p,'sym')   % DJ, check class type, 12-10-2021
    % If the upper and lower bounds are close (absolute or relative), return them
        ach = double(subs(p,num2cell(vars).',num2cell(xopt).'));
        
        if min(abs(ach/GAM-1),abs(ach-GAM)) > 1e-4 ;
            xopt = [];
            return;
        end
        
        % Check inequality and equality constraints
        if length(ineq)>1
            ineq = double(subs(ineq,num2cell(vars).',num2cell(xopt).'));
        else
            ineq = 0;
        end;
        if length(eq)>1
            eq = double(subs(eq,num2cell(vars).',num2cell(xopt).'));
        else
            eq = 0;
        end;
        if min(ineq)<-1e-6 | max(abs(eq))>1e-6
            xopt = [];
            return;
        end;
    elseif isa(p,'polynomial')
    % If the upper and lower bounds are close (absolute or relative), return them
        ach = double(subs(p,vars,xopt));
        
        if min(abs(ach/GAM-1),abs(ach-GAM)) > 1e-4
            xopt = [];
            return;
        end
        
        % Check inequality and equality constraints
        if length(ineq)>1
            ineq = double(subs(ineq,vars,xopt));
        else
            ineq = 0;
        end
        if length(eq)>1
            eq = double(subs(eq,vars,xopt));
        else
            eq = 0;
        end
        if min(ineq)<-1e-6 || max(abs(eq))>1e-6
            xopt = [];
            return;
        end        
    else
    % If the upper and lower bounds are close (absolute or relative), return them
        psol = sosgetsol(prog,p,16);
        ach = double(subs(psol,vars,xopt));
        
        if min(abs(ach/GAM-1),abs(ach-GAM)) > 1e-4
            xopt = [];
            return;
        end
        
        % Check inequality and equality constraints
        if length(ineq)>1
            ineq = double(subs(ineq,vars,xopt));
        else
            ineq = 0;
        end
        if length(eq)>1
            eq = double(subs(eq,vars,xopt));
        else
            eq = 0;
        end
        if min(ineq)<-1e-6 || max(abs(eq))>1e-6
            xopt = [];
            return;
        end        
    end
    
end;