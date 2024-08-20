function [info,dopt,sossol] = gsosopt(pconstr,x,t,opts)
% function [info,dopt,sossol] = gsosopt(pconstr,x,t,opts)
%
% DESCRIPTION
%   This function solves the following bilinear optimization problem with
%   polynomial equality and sum-of-squares constraints.
%        min_{d,t} t
%        subject to:
%           gi(x,d) + t*fi(x,d) is sum of squares (i=1,...,Ns)
%           fi(x,d) is sum of squares (i=1,...,Ns)
%           pj(x,d)==0                (j=1,...,Ne)
%   where all polynomials are affine in the decision variables d. 
%   If fi(x,d)=0 then the i^th constraint is a linear SOS constraint but 
%   if fi(x,d) depends on d then the constraint is bilinear in t and d.  
%   This optimization is quasiconvex and is solved with bisection.
%
% INPUTS
%   pconstr: Np-by-1 array of equality and SOS constraints. Each constraint
%       must have the form gi(x,d)+t*pi(x,d) where gi and pi are affine in 
%       the decision variables d. For each SOS constraint where pi(x,d) is 
%       nonzero the user must specify another constraint that enforces 
%       pi(x,d) in SOS.
%   x: Nx-by-1 array of pvars. These are the polynomial variables in the
%       optimization.  All remaining variables are decision variables.
%   t: 1-by-1 pvar. t is the objective function in the optimization.
%   opts: Options for optimization.  See GSOSOPTIONS for more details.
%
% OUTPUTS
%   info: Structure with fields feas, tbnds, opts, sosdata, sos2sdp, 
%       sdpdata, sdpsol, and sosfeasfh. feas=1 if the SOS optimization is 
%       feasible for some t in [minobj,maxobj] and feas=0 otherwise. tbnds 
%       is a 1-by-2 vector [tlb, tub] giving a lower bound tlb and upper 
%       bound tub on the minimum value of t. tbnds is empty if feas=0.
%       The remaining fields contain data about the billinear SOS 
%       optimization and its related SDP.
%   dopt: (Nd+1)-by-2 polynomial array of the optimal decision variables
%       (d,t). dopt(:,1) contains the decision variables and dopt(:,2) 
%       contains the optimal values. Use the subs command to replace the
%       decision variables in any polynomial with their optimal values,
%       e.g. the minimal value of t can be computed by subs(t,dopt).
%       dopt is returned as if empty if the optimization is infeasible.
%   sossol: Np-by-1 structure array with fields p, z, and Q. sossol(i).p
%       is pconstr(i) evaluated at the optimal decision variables.
%       sossol(i).z and sossol(i).Q are the vector of monomials and
%       positive semidefinite matrix for the Gram matrix decomposition of
%       sossol(i).p. sossol is empty if the optimization is infeasible.
%
% SYNTAX
%  [info,dopt,sossol] = gsosopt(pconstr,x,t)
%  [info,dopt,sossol] = gsosopt(pconstr,x,t,opts)
%    Solve a Generalized SOS minimization with or without options
%
% See also sosopt, gsosoptions, polydecvar, sosdecvar


% 12/04/09  PJS     Initial Coding
% 11/01/10  PJS     Update for polyconstr, sosoptions, == constraints


% XXX Improve help: Document syntax w/ non-scalar opts.minobj and opts.maxobj

%---------------------------------------------------------------------
% Error checking
%---------------------------------------------------------------------

% Check that x is a pvar and convert pvars to a cell array of strings
if ispvar(x)
    x = x.var;
elseif ~iscellstr(x)
    error('x must be a vector of pvars or a cell array of strings');
end

% Check that t is a pvar
if isa(t,'char') && any(size(t)==1)
    t = pvar(t);
elseif ~(ispvar(t) && length(t)==1)
    error('t must be a pvar or string.');
end
tchar = t.varname{1};

%---------------------------------------------------------------------
% Process Inputs
%---------------------------------------------------------------------

% Get options
if nargin <= 3
    opts = gsosoptions;
end
if strcmp(opts.checkfeas,'off')
    warning(['GSOSOPT bisection is based on a feasibility check. ' ...
        'opts.checkfeas = ''off'' is not allowed. ' ...
        'Setting opts.checkfeas to ''fast''.']);
    opts.checkfeas = 'fast';
end
if ~strcmp(opts.form,'image')
    warning(['Currently only the image form implemented for GSOSOPT. ' ...
        'Setting opts.form to ''image''.']);
    opts.form = 'image';
end

tmin = opts.minobj;
tmax = opts.maxobj;
absbistol = opts.absbistol;
relbistol = opts.relbistol;
    
% Convert polynomial/cell pconstr (old syntax) to polyconstr object
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
        tmp(i1) = polynomial(pconstr{i1}) >= 0;
    end
    pconstr = tmp(:);
elseif ~isa(pconstr,'polyconstr')
    error('Invalid syntax: first input argument must be a polyconstr object');
end

% Find decision variables
allvar = [];
for i1=1:Np
    allvar = [allvar; pconstr(i1).LeftSide.varname];
    allvar = [allvar; pconstr(i1).RightSide.varname];
end
if iscellstr(x)
    Nx = length(x);
    x = polynomial(speye(Nx),speye(Nx),x,[Nx 1]);
else
    x = polynomial(x);
end
%decvars = setdiff(allvar,x.var);  % setdiff returns a unique list

%------------------------------------------------------------------
% Pre-processing for fast SDP set-up
%------------------------------------------------------------------

% XXX Check here that constraints are linear in decvars?
% (Otherwise SOSOPT might error out and give errors with quasilinear vars)

% Create quasi-linear constraints in vars (d,dbar) where dbar:= t*d
quasilinc = polyconstr;
dbarprefix = 'tGSOSOPTINTERNAL';
zerop = polynomial(0);
for i1=1:Np
    % Convert to one-sided constraint, s(x,d)==0 or s(x,d)>=0
    RelOp = pconstr(i1).RelOp;
    s = pconstr(i1).OneSide;
    
    % Split s(x,d,t) into the form g(x,d)+t*p(x,d)
    [g,p,h] = collect(  s , setdiff(s.varname,t.varname) );
    if ~( isempty(h) || isequal(h,t) ) 
        error(['Constraints must be affine in ' tchar]);
    end
    if isequal(g,zerop) && ~strcmp(RelOp,'==')
        % If g=0 then the constraint is t*p>=0. This is either
        % redundant or in conflict with the default constraint p>=0.
        % XXX Improve error message
        error('g(x,d) must be nonzero in each constraint');
    end
    
    % Split: p(x,d) = a(x) + b(x)*d
    % Define dbar=t*d and then t*p(x,d) = t*a(x)+b(x)*dbar
    %  -->s(x,d)=g(x,d)+p(x,dbar) is affine in (d,dbar)
    if isempty(p)
        pbar = 0;
    else
        [a,b,dbar]=collect(p,x);
        pbar = pvar(dbarprefix)*a;
        if ~isempty(dbar)
            dbarvar = dbar.varname;
            dbarvar = [repmat([dbarprefix '_'],[length(dbarvar) 1]) char(dbarvar)];
            dbar.varname = cellstr(dbarvar);
            %pbar = pvar(dbarprefix)*a+b*dbar;
            pbar = pbar + b*dbar;
        end
    end
    
    if strcmp(RelOp,'==')
        quasilinc = [quasilinc; g+pbar==0];
    else
        quasilinc = [quasilinc; g+pbar>=0];        
        % XXX Automatically enforce p in SOS?
        % quasilinc = [quasilinc; p>=0];
    end
end

% Check that user did not have a variable that conflicts with dbarprefix
if ~isempty( strmatch(dbarprefix,allvar) )
    errstr = [dbarprefix ' is used internally by gsosopt. No' ...
        ' variable names may start with this prefix.'];
    error(errstr);
end

% Use SOSOPT to convert the quasi-linear SOS problem into Sedumi form.
% This will be used for a fast setup at each bisection step.
obj = [];
setupopts = opts; % XXXX
setupopts.solver = 'setup';
setupopts.simplify = 'off';   
info=sosopt(quasilinc,x,obj,setupopts);
info.opts = opts;
info.sosdata.obj = t;
info.sosfeasfh = @(info,ttry) LOCALsosfeas(info,ttry,dbarprefix);
info.tbnds = [];
info.iter  = [];
info = rmfield(info,'obj');


if strcmp(opts.display,'on') || strcmp(opts.display,'pcontain')
    fprintf('\n----------- Beginning Bisection\n');
end

%------------------------------------------------------------------
% Initialize upper bound from tmax
%------------------------------------------------------------------
dopt = [];
sossol = [];
tlb = -inf;
tub = inf;
Ntmax = length(tmax);
if Ntmax == 1
    tub = tmax;
else
    % Undocumented: If tmax is a vector then the code below will try to
    % check all the values listed in tmax until it finds a feasible
    % solution.  This can be used if user wants to speed up the bisection
    % search with a good initial guess for tub.
    tmax = sort(tmax,'ascend');
    for i1 = 1:Ntmax
        ttry = tmax(i1);
        if nargout==3
            [tinfo,d,sol] = LOCALsosfeas(info,ttry,dbarprefix);
        else
            [tinfo,d] = LOCALsosfeas(info,ttry,dbarprefix);
            sol = [];
        end
        feas = tinfo.feas;
        
        if strcmp(opts.display,'on');
            fprintf(['Verify Upper Bound: ' tchar '_try = %4.3f\t'],ttry)
            fprintf('Feas = %d\n',feas);
        elseif strcmp(opts.display,'pcontain');
            fprintf(['Verify Lower Bound: ' tchar '_try = %4.3f\t'],-ttry)
            fprintf('Feas = %d\n',feas);
        end
        if feas==1
            tub = ttry;
            dopt = d;
            sossol = sol;
            info.sdpsol = tinfo.sdpsol;
            break;
        else
            tlb = ttry;
        end
    end
    if feas==0
        % Bisection was invalid even for largest val of tmax
        info.tbnds = [];
        dopt = [];
        sossol = [];
        return;
    end
end

%------------------------------------------------------------------
% Intialize lower bound from tmin
%------------------------------------------------------------------
tmin( tmin<=tlb | tmin>=tub ) = [];
if ~isinf(tlb)
    tmin = [tmin(:); tlb];
end
Ntmin = length(tmin);
if Ntmin ==1
    tlb = tmin;
elseif ~isequal(tmin,tlb)
    % Undocumented: If tmin is a vector then the code below will try to
    % check all the values listed in tmin until it finds an infeasible
    % solution.  This can be used if the user wants to speed up the
    % bisection search with a good initial guess at tmin.
    tmin = sort(tmin,'descend');
    for i1 = 1:Ntmin
        ttry = tmin(i1);
        if nargout==3
            [tinfo,d,sol] = LOCALsosfeas(info,ttry,dbarprefix);
        else
            [tinfo,d] = LOCALsosfeas(info,ttry,dbarprefix);
            sol = [];
        end
        feas = tinfo.feas;
        
        if strcmp(opts.display,'on');
            fprintf(['Verify Lower Bound: ' tchar '_try = %4.3f\t'],ttry)
            fprintf('Feas = %d\n',feas);
        elseif strcmp(opts.display,'pcontain');            
            fprintf(['Verify Upper Bound: ' tchar '_try = %4.3f\t'],-ttry)
            fprintf('Feas = %d\n',feas);            
        end
        if feas==1
            tub = ttry;
            dopt = d;
            sossol = sol;
            info.sdpsol = tinfo.sdpsol;
        else
            tlb = ttry;
            break;
        end
    end
    if feas==1
        % Bisection was valid even for smallest val of tmin
        info.tbnds = [ttry ttry];
        return;
    end
end

%------------------------------------------------------------------
% Perform Bisection
%------------------------------------------------------------------
info.iter = 0;  % Ensure at least one pass through
info.sdpsol = [];
%while (tub-tlb>absbistol && tub-tlb>relbistol*tlb) || (go==1)
while (tub-tlb>absbistol && tub-tlb>relbistol*abs(tlb)) || ~info.iter
    % Set gamma level for feasibility problem
    ttry = (tub+tlb)/2;
    
    % Solve feasibility problem
    if nargout==3
        % XXX Why form sol at each step? Do it only at the end?
        % We need to do it if checkfeas = 'full'
        [tinfo,d,sol] = LOCALsosfeas(info,ttry,dbarprefix);
    else
        [tinfo,d] = LOCALsosfeas(info,ttry,dbarprefix);
        sol = [];
    end
    feas = tinfo.feas;
    %info.sdpsol = tinfo.sdpsol;
    
    % Display Bisection Status
    if strcmp(opts.display,'on');
        fprintf([tchar 'lb  = %4.5f \t ' tchar 'try = %4.5f \t'],tlb,ttry);
        fprintf([tchar 'ub = %4.5f \t Feas = %d \n'],tub,feas);
    elseif strcmp(opts.display,'pcontain');
        fprintf([tchar 'lb  = %4.5f \t ' tchar 'try = %4.5f \t'],-tub,-ttry);
        fprintf([tchar 'ub = %4.5f \t Feas = %d \n'],-tlb,feas);            
    end
    info.sdpsol = [info.sdpsol tinfo.sdpsol];
    % Update bisection bounds
    if feas==1
        tub = ttry;
        dopt = d;
        sossol = sol;
      
    else
        tlb = ttry;
    end
    
    info.iter = info.iter+1;
end

% Set output variable with upper/lower bounds
if isempty(dopt)
    info.feas = 0;
    info.tbnds = [];
else
    info.feas = 1;
    info.tbnds = [tlb tub];    
end
info = orderfields(info,{'feas','tbnds','opts','sosdata','sos2sdp',...
    'sdpdata','sdpsol','sosfeasfh','iter'});

return;

%-----------------------------------------------------------
% Local Function:
%   SOS feasibility problem: Given t, find d such that
%    si(x,d) + t*pi(x,d) is sum of squares (i=1,...,Ns)
%    pi(x,d) is sum of squares (i=1,...,Ns)
%    A*d<=b
%-----------------------------------------------------------
function [tinfo,dopt,sossol] = LOCALsosfeas(info,ttry,dbarprefix)

% Get SDP info
opts = info.opts;
decvars = info.sosdata.decvars;
dv2x  = info.sos2sdp.dv2x;
t = info.sosdata.obj;
A = info.sdpdata.A;
b = info.sdpdata.b;
K = info.sdpdata.K;

% Find indices of dbar vars and corresponding true dec vars
L = length(dbarprefix);
dbaridx = strmatch(dbarprefix,decvars);
dbar = decvars(dbaridx);
d = char(dbar); 
d = cellstr(d(:,L+2:end));
[tf1,idx1]=ismember(d,decvars);
[tf2,idx2]=ismember(dbar,decvars);

% Handle vars that enter bilinearly as t*d 
blidx = find(idx1==0);
tmpidx = idx2(blidx);
A(:,dv2x(tmpidx)) = ttry*A(:,dv2x(tmpidx));
decvars(dbaridx(blidx)) = d(blidx);
dbar(blidx) = [];
dbaridx(blidx) = [];
idx1(blidx) = [];
idx2(blidx) = [];
Ndbar = length(dbar);

% Need to handle terms linear in t (as opposed to bilinear in t & d)
%   idx1: Locations of d in decvars
%   idx2: Locations of dbar = t*d in decvars
%   idx3: Location of t in decvars
%   dbaridx: Locations of (t,dbar) in decvars
tmpidx = strcmp(dbarprefix,dbar);
%d(tmpidx) = [];
idx1(tmpidx) = [];
idx3 = idx2(tmpidx);
idx2(tmpidx) = [];

% Create linear SOS problem at the current value of ttry
A(:,dv2x(idx1)) = A(:,dv2x(idx1))+ttry*A(:,dv2x(idx2));
if ~isempty(idx3)
    b = b - ttry*A(:,dv2x(idx3));
end
A(:,dv2x(dbaridx)) = [];
c = sparse(size(A,2),1);
K.f = K.f - Ndbar;

% Update decvar list and dv2x indexing
decvars(dbaridx) = [];
Ndv = length(decvars);

tmpidx = 1:size(A,2);
tmpidx(dv2x(dbaridx))=[];
dv2x(dbaridx) = [];

xold2xnew(tmpidx) = 1:length(tmpidx);
dv2x = xold2xnew(dv2x);

% Do SOS simplification
z = info.sos2sdp.z;
Ncsosdv = length( info.sosdata.standform.sosdecvar.idx );
if strcmpi(opts.simplify,'on');
    [A,b,K,z,dv2x,Nfv,feas,zrem] = sospsimplify(A,b,K,z,dv2x,Ncsosdv);
else
    zrem = cell(length(K.s),1);
    feas = 1;
end
sdpdata2.A = A;
sdpdata2.b = b;
sdpdata2.c = c;
sdpdata2.K = K;

% Store info associated with this value of t
tinfo = rmfield(info,{'sosfeasfh','tbnds'});
tinfo.feas =feas;
tinfo.t = ttry;
tinfo.sdpdata = sdpdata2;

tinfo.sos2sdp.z = z;
tinfo.sos2sdp.zremoved = zrem;
tinfo.sos2sdp.dv2x = dv2x;
% tinfo.sos2sdp.Nfv = Nfv;
tinfo = orderfields(tinfo,{'feas','t','opts','sosdata','sos2sdp','sdpdata','sdpsol','iter'});

% Exit if sossimplify detected infeasibility
if feas
    if ~isempty(A)
        [xs,ys,solverinfo] = solvesossdp(sdpdata2,opts);
    else
        % 6/22/10: Abhijit found an example that was feasible but
        % SOSSIMPLIFY removed all variables.
        xs = zeros(size(A,2),1);
        ys = zeros(size(A,1),1);
        solverinfo.iter = 0;
        solverinfo.feasratio = 1;
        solverinfo.pinf = 0;
        solverinfo.dinf = 0;
        solverinfo.numerr = 0;
        solverinfo.timing = 0;
        solverinfo.cpusec = 0;
    end
else
    % Problem was determined to be infeasible by SOSSIMPLIFY
    xs = zeros(size(A,2),1);
    ys = zeros(size(A,1),1);
    solverinfo.iter = 0;
    solverinfo.feasratio = -1;
    solverinfo.pinf = 1;
    solverinfo.dinf = 1;
    solverinfo.numerr = 0;
    solverinfo.timing = 0;
    solverinfo.cpusec = 0;
end

% sdpsol: x,y, solverinfo
sdpsol.x = xs;
sdpsol.y = ys;
sdpsol.solverinfo = solverinfo;
tinfo.sdpsol = sdpsol;

%---------------------------------------------------------------------
% Create sossol and dopt for output
%---------------------------------------------------------------------

% Create dopt
if Ndv ~= 0
    dvals = zeros(Ndv,1);
    idx = find(dv2x~=0);
    dvals( idx ) = xs( dv2x(idx) );
    %dopt = [decvars(:) num2cell(dvals)];
    %dopt = [char(t), ttry; dopt];
    dopt = [polynomial(decvars(:)) dvals];
    dopt = [t, ttry; dopt];
        
    dbaropt = polynomial(zeros(Ndv,2));
    for i1=1:Ndv
        di = char(dopt(i1+1,1));
        dbaropt(i1,1) = polynomial({[dbarprefix '_' di{1}]});
        dbaropt(i1,2) = ttry*dopt(i1+1,2);
    end
    dbaropt(Ndv+1,1) = polynomial({dbarprefix});
    dbaropt(Ndv+1,2) = ttry;
else
    dopt = [];
end

% Create sossol
sossol = [];
if nargout==3 || ~strcmp(opts.checkfeas,'fast')
    cnt = Nfv;
    dall =[dopt;dbaropt];
    sosdata = info.sosdata;
    pconstr = sosdata.pconstr;
    Np = length(pconstr);
    
    scaling = sosdata.scaling;
    standform = sosdata.standform;
    Nceq = length(standform.equality.idx);
    Nclp = length(standform.linprog.idx);
    Ncsosdv = length(standform.sosdecvar.idx);
    Ncsos = length(standform.sos.idx);
    
    for i1=1:Np
        % Evaluate pconstr at dopt
        p = pconstr(i1);
        if Ndv~=0
            p = subs(p,dall);
        end
        sossol(i1,1).p = p;
    end
    
    % Fill sossol for equality constraints
    % XXX Return anything for equality constraints?
    idx=standform.equality.idx;
    for i1=1:Nceq
        sossol(idx(i1)).z = [];
        sossol(idx(i1)).Q = [];
    end

    % Fill sossol for LP constraints
    idx = standform.linprog.idx;
    if Nclp>0
        sz = length(standform.linprog.data{2});
        xidx = cnt+(1:sz);
        Qvec = xs(xidx);
        for i1=1:Nclp
            % Unscale: XXX Currently no LP scaling.
            sossol(idx(i1)).z = 1;
            sossol(idx(i1)).Q = Qvec(i1);
        end
        if ~isempty(xidx)
            cnt = xidx(end);
        end
    end
    
    % Fill sossol for SOS constraints
    for i1=1:Ncsosdv+Ncsos
        if i1<=Ncsosdv
            idx=standform.sosdecvar.idx(i1);
            scl = 1;
        else
            idx=standform.sos.idx(i1-Ncsosdv);
            scl = scaling.sos(i1-Ncsosdv);
        end
        sz = K.s(i1);
        xidx = cnt+(1:sz^2);
        Q = reshape( xs(xidx), sz, sz );
        
        % Unscale solution
        if scl~=1
            Q = Q*scl;
        end
        
        % Store Gram Matrix solution for output
        sossol(idx).z = z{i1};
        sossol(idx).Q = Q;
        
        if ~isempty(xidx)
            cnt = xidx(end);
        end
    end    
end

% Feasibility Check
feas = checkfeas(tinfo,sossol,opts);
tinfo.feas = feas;

