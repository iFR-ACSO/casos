function  [info,dopt,sossol] = sosopt(pconstr,x,obj,opts)
% function [info,dopt,sossol] = sosopt(pconstr,x,obj,opts)
%
% DESCRIPTION
%   This function solves an optimization or feasibility problem with
%   polynomial equality and sum-of-squares constraints.
%           min c'*d
%           subject to:   si(x,d) is sum of squares  (i=1,...,Ns)
%                         pj(x,d) == 0               (j=1,...,Ne)
%   d is a vector of decision variables. Refer to the polydecvar and
%   sosdecvar functions for the creation of polynomial decision variables.
%
% INPUTS
%   pconstr: Np-by-1 array of equality and SOS constraints. All
%       polynomials must be linear functions of the decision variables.
%   x: Nx-by-1 array of pvars. These are the polynomial variables in the
%       optimization.  All remaining variables are decision variables.
%   obj: 1-by-1 polynomial. obj is the objective function and must be a
%       linear function of the decision variables. Set obj =[] for a
%       feasibility problem
%   opts: Options for optimization.  See SOSOPTIONS for more details.
%
% OUTPUTS
%   info: Structure with fields feas, opts, sosdata, sos2sdp, sdpdata, and
%       sdpsol. feas=1 if the SOS optimization is feasible and feas=0
%       otherwise. obj is the minimal value of the objective function.
%       obj is +inf if feas=0. The remaining fields contain data about 
%       the SOS optimization and the related SDP.
%   dopt: Nd-by-2 polynomial array of the optimal decision variables. 
%       dopt(:,1) contains the decision variables and dopt(:,2) contains
%       the optimal values. Use the subs command to replace the decision
%       variables in any polynomial with their optimal values, e.g. 
%       subs(obj,dopt) returns the minimal value of the objective function.
%       dopt is returned as if empty if the optimization is infeasible.
%   sossol: Np-by-1 structure array with fields p, z, and Q. sossol(i).p
%       is pconstr(i) evaluated at the optimal decision variables.
%       sossol(i).z and sossol(i).Q are the vector of monomials and
%       positive semidefinite matrix for the Gram matrix decomposition of
%       sossol(i).p. sossol is empty if the optimization is infeasible.
%
% SYNTAX
%  [info,dopt,sossol] = sosopt(p)
%    Determine if p is SOS. See also ISSOS.
%  [info,dopt,sossol] = sosopt(pconstr,x);
%  [info,dopt,sossol] = sosopt(pconstr,x,opts);
%    Solve a feasibility problem with or without options
%  [info,dopt,sossol] = sosopt(pconstr,x,obj);
%  [info,dopt,sossol] = sosopt(pconstr,x,obj,opts);
%    Solve a SOS minimization with or without options
%
% See also issos, sosoptions, polydecvar, sosdecvar

% 3/05/08  PJS     Initial Coding
% 4/25/09  PJS     Cleaning up documentation + SOS dec var handling
% 5/6/09   PJS     Modify to allow set-up call
% 10/28/10 PJS     Update for polyconstr, sosoptions, == constraints

% XXX Should I return dopt as an Nd-by-2 poly array with the
% first column containing the dec vars (as pvars) and the second
% column containing their values?

%---------------------------------------------------------------------
% Error checking
%---------------------------------------------------------------------
error(nargchk(1,4,nargin))

%---------------------------------------------------------------------
% Process Inputs
%---------------------------------------------------------------------
if nargin < 2
    x = []; obj = []; opts = sosoptions;
elseif nargin < 3
    if isa(x,'sosoptions')
        opts = x;
        x = []; obj = [];
    else
        obj = []; opts = sosoptions;
    end
elseif nargin <4
    if isa(obj,'sosoptions')
        opts= obj;
        obj = [];
    else
        opts = sosoptions;
    end
end
if ~isa(opts,'sosoptions')
    error('Options must be an sosoptions object');
end

% Convert polynomial/cell pconstr (old syntax) to polyconstr object
x = polynomial(x);
pconstr = pconstr(:);
Np = length(pconstr);
if isa(pconstr,'double') || isa(pconstr,'polynomial')
    pconstr = polynomial(pconstr);
    tmp = polyconstr;
    for i1=1:Np
        tmp(i1) = pconstr(i1) >= 0;
    end
    pconstr = tmp(:);
elseif iscell(pconstr)
    tmp = polyconstr;
    for i1=1:Np
        if isa(pconstr{i1},'polyconstr')
            tmp(i1) = pconstr{i1};
        elseif isa(pconstr{i1},'polynomial')
            % Convert to a polyconstr
            tmp(i1) = polynomial(pconstr{i1}) >= 0;
        else
            error('If pconstr is a cell array then contents must be a polynomial or polyconstr object');
        end
    end
    pconstr = tmp(:);
elseif ~isa(pconstr,'polyconstr')
    error('Invalid syntax: first input argument must be a polyconstr object');
end

%---------------------------------------------------------------------
% Find decision variables
%---------------------------------------------------------------------
allvar = [];
for i1=1:Np
    allvar = [allvar; pconstr(i1).LeftSide.varname];
    allvar = [allvar; pconstr(i1).RightSide.varname];
end

if isempty(x)
    Nx = length(allvar);
    x = polynomial(speye(Nx),speye(Nx),allvar,[Nx 1]);
elseif iscellstr(x)
    Nx = length(x);
    x = polynomial(speye(Nx),speye(Nx),x,[Nx 1]);
end
decvars = setdiff(allvar,x.var);  % setdiff returns a unique list

%---------------------------------------------------------------------
% Map pconstr constraints to a standard form
%
% Standard form types are 'equality', 'linprog', 'sosdecvar', or 'sos'.
% This classification is used for efficient optimization implementation.
%---------------------------------------------------------------------
standform = [];

% Convert to one-sided constraint, s(x,d)==0 or s(x,d)>=0
s = pconstr.OneSide;

% equality: s(x,d) == 0
eqidx = find( strcmp('==',pconstr.RelOp) );
Neq = length(eqidx);
standform.equality.idx = eqidx;
standform.equality.data = s(eqidx);

% Linear Progamming constraints
% Note: Treat constraints of the form d>=0 as sosdecvar constraints
ineqidx = ~strcmp('==',pconstr.RelOp);
[isdec,xidx]=ismember(x.var,s.var);
xidx( xidx==0 ) = [];
xtermidx = find( sum( s.degmat(:,xidx) , 2 )~=0 );
lpidx = sum( abs(s.coef(xtermidx,:)) , 1)==0;
lpidx = find( lpidx(:)  & ineqidx );
[b,A,d] = collect(s(lpidx),[]);
b = double(b);
A = -double(A);

if ~isempty(b)
    % linprog: s(d) = b-A*d>=0
    idx = b==0 & sum( A~=0 ,2)==1 & sum( A ,2)==-1;
    standform.linprog.idx = lpidx(~idx);
    standform.linprog.data = {A(~idx,:),b(~idx),d};

    % sosdecvar: s(x,d) = z(x)'*D*z(x)>=0
    idx = find(idx);
    standform.sosdecvar.idx = lpidx(idx);
    standform.sosdecvar.data = cell(length(idx),3);
    for i1=1:length(idx)
        D = s(lpidx(idx(i1)));
        standform.sosdecvar.data(i1,1:3) =  {1,D,D};  %  = {z,D,s}
    end
else
    standform.linprog.idx = [];
    standform.linprog.data = cell(0,3);
    standform.sosdecvar.idx = [];
    standform.sosdecvar.data = cell(0,3);
end
    
% SOS Constraints
sosidx = setdiff(1:Np,[eqidx;lpidx]);
standform.sos.idx = [];
for i1=1:length(sosidx)
    idx = sosidx(i1);
    [D,z] = detectsosdecvar(s(idx),x);
    if ~isempty(D)
        % sosdecvar: s(x,d) = z(x)'*D*z(x)>=0
        standform.sosdecvar.idx = [standform.sosdecvar.idx; idx];
        standform.sosdecvar.data = [standform.sosdecvar.data; {z,D,s(idx)}];
    else
        % sos: s(x,d) >=0
        standform.sos.idx = [standform.sos.idx; idx];
    end
end
standform.sos.data = s( standform.sos.idx );
Nsos = length(standform.sos.idx);

%---------------------------------------------------------------------
% Scale (Equality and SOS) Constraints
% XXX Should we scale LP constraints?
%---------------------------------------------------------------------
eqscl = ones( Neq ,1 );
sosscl = ones( Nsos ,1 );
if strcmpi(opts.scaling,'on');
    for i1=1:Neq
        s = standform.equality.data(i1);
        eqscl(i1) = sqrt(norm(s.coefficient));
        s.coefficient = s.coefficient/eqscl(i1);
        standform.equality.data(i1) = s;
    end
    for i1=1:Nsos
        s = standform.sos.data(i1);
        sosscl(i1) = sqrt(norm(s.coefficient));
        s.coefficient = s.coefficient/sosscl(i1);
        standform.sos.data(i1) = s;
    end    
end
scaling.equality = eqscl;
scaling.sos = sosscl;

%---------------------------------------------------------------------
% Store all data about the SOS problem formulation
%---------------------------------------------------------------------
sosdata.pconstr = pconstr;
sosdata.standform = standform;

sosdata.x = x;
sosdata.decvars = decvars;
sosdata.obj = obj;
sosdata.scaling = scaling;

info.sosdata = sosdata;
info.opts = opts;

%---------------------------------------------------------------------
% Call Image/Kernel Forms of SOS Optimization
%---------------------------------------------------------------------
if (nargout==3) || any( strcmpi(opts.checkfeas,{'full';'both'}) )
    if strcmpi(opts.form,'image');
        [info,dopt,sossol]=sosp(sosdata,opts);
    elseif strcmpi(opts.form,'kernel');
        [info,dopt,sossol]=sosd(sosdata,opts);
    end
else
    if strcmpi(opts.form,'image');
        [info,dopt]=sosp(sosdata,opts);
    elseif strcmpi(opts.form,'kernel');
        [info,dopt]=sosd(sosdata,opts);
    end
    sossol = [];
end

%---------------------------------------------------------------------
% Check feasibility of solution
%---------------------------------------------------------------------
feas = checkfeas(info,sossol,opts);
info.feas = feas;
info = orderfields(info,{'feas','obj','opts','sosdata','sos2sdp','sdpdata','sdpsol'});
if feas==0
    % Infeasible solution
    info.obj = +inf;
    sossol = [];
    dopt = [];
end

