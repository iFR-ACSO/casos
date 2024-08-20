function  [info,dopt,sossol] = sosp(sosdata,opts)
% function [info,dopt,sossol] = sosp(sosdata,opts)
%
% DESCRIPTION
%   This function solves the SOS optimization or feasibility problem in
%   image (primal) form.  An SOS optimization contains constraints of the
%   form (see sosopt for more details):
%               s(x,d) is sum of squares
%               p(x,d) = 0
%   where d are the decision variables and x are the polynomial variables.
%
%   This function solves the optimization by introducing a matrix decision
%   variable Q for each SOS constraint and rewriting the constraint as:
%               s(x,d)=z(x)'*Q*z(x)
%               Q>=0
%   The constraint s(x,d)=z(x)'*Q*z(x) is then represented as a linear
%   equality constraint between the entries of Q, the decision variables
%   d, and the numeric coefficients of s. The polynomial equality 
%   constraints are converted into equality constraints on the decision
%   variables.
%
% See also sosopt, sosoptions, sosd

% 1/29/08  PJS     Initial Coding
% 2/26/08  PJS     Added interface to DSDP
% 11/11/08 PJS     Separately handle pure sos variables
% 4/25/09 PJS      Clean up documentation and code for pure sos vars
% 5/6/09 PJS       Modify to allow set-up call
% 5/29/09 PJS      Modify to handle symmetric SOS coef matrix
% 9/21/09 PJS      Modify for monomial and decision variable reduction
% 10/28/10 PJS     Update for polyconstr, sosoptions, == constraints


%---------------------------------------------------------------------
% Formulate SOS Optimization as an SDP
%---------------------------------------------------------------------

% Get data about SOS problem formulation
pconstr = sosdata.pconstr;
Np = length(pconstr);

standform = sosdata.standform;
Nceq = length(standform.equality.idx);
Nclp = length(standform.linprog.idx);
Ncsosdv = length(standform.sosdecvar.idx);
Ncsos = length(standform.sos.idx);

x = sosdata.x;
decvars = sosdata.decvars;
obj = sosdata.obj;
scaling = sosdata.scaling;

% Compute indexing, dv2xs, of decision vars into sedumi vars, xs
% (see definition of xs below)
Nlpv = Nclp;
dv2xs = [];
sosvars = [];
Nsosv = 0;
for i1=1:Ncsosdv
    % D is symmetric.  Grab dec vars from lower half of D
    % XXX Revisit: What if the same sos poly is listed twice in
    % the constraint list?  decvars is a unique list of vars but
    % it will be repeated twice here.
    cdata = standform.sosdecvar.data(i1,:);
    D = cdata{2};
    lz = size(D,1);
    idx = find(tril(ones(lz)));
    sosvars = [sosvars; char(D(idx))];
        
    % Index mapping from symmetric decvars to non-symmetric
    dv2xs = [dv2xs Nsosv+idx(:)'];
    Nsosv = Nsosv+lz^2;
end
freevars = setdiff(decvars,sosvars);
Nfv = length(freevars);
dv2xs = [1:Nfv  dv2xs+Nfv+Nlpv];
decvars = [freevars; sosvars];
Ndv = length(decvars);

%---------------------------------------------------------------------
% Define Objective and SOS/LP Constraints in Sedumi Primal Form
%       min cs'*xs   s.t. As*xs = bs  ;  xs >= 0
%
% Equality Constraints: The constraint "s(x,d)==0" can be expressed
%   as Ad*d = b for appropriately defined vector b and matrix Ad.
% SOS Constraints: The constraint "s(x,d) is SOS" is equivalent to the
%   constraints s(x,d) = z(x)'*Q*z(x) and Q>=0. This can be expressed
%   as Ad*d + Aq*Q(:) = b and Q>=0 for appropriately defined vector b
%   and matrices Ad and Aq. Q is a PSD slack matrix.
% LP Constraints: A*d<=b is defined as A*d+y=b and y>=0 where
%   y is a slack variable.
% SOS Variables: As described above, generic SOS constraints are
%   represented as a set of equality constraints with a slack PSD
%   matrix.  SOS variables do not require the introduction of a slack
%   matrix.  They are directly reresented as z(x)'*M*z(x) where M is a
%   non-symmetric matrix of variables [due to Sedumi format) and z(x)
%   is a vector of monomials. This constraint can be directly represented
%   as M>=0. SOS variables are represented using a symmetric matrix of
%   decision variables. These symmetric matrix decision vars map into the
%   lower half of the Sedumi nonsymmetric matrix vars.
%
% The order of variables in the sedumi format is free variables, LP
% variables, and semidefinite variables. For the SOS program, this
% ordering is:
%   xs = [Free decision vars (Nfv-by-1);
%         y = LP slack vars  (Nlpv-by-1);
%         M(:) = SOS dec vars in the form of PSD matrices (Nsosv-by-1);
%         Q(:) = PSD slack matrices for SOS inequalities]
%---------------------------------------------------------------------

% Define Number of Free Decision Variables
K.f = Nfv;

% Define Equality constraints: A*d=b
As = sparse(0,Nfv+Nlpv+Nsosv);
bs = sparse(0,1);
for i1=1:Nceq
    % Get equality constraint data    
    s = standform.equality.data(i1);
    
    % Collect s(x,d) into the form g0(x)+g1(x)*d1 + ... + gN(x)*dN
    % where di are decision variables.
    [g0,g,d] = collect(s,x);
    if ispvar(d)
        % Find the decision variables involved in the equality constraint
        [isdec,didx]=ismember(char(d),decvars);
    elseif isempty(d)
        isdec=0;
        didx=[];
    else
        % Error: Nonlinear dependence on decision variables
        tmpidx = find( sum(d.degmat,2) > 1);
        dterm = polynomial(1,d.degmat(tmpidx(1),:),d.varname,[1 1]);
        error( LOCALerrstr(dterm,standform.equality.idx(i1)) )
    end
    
    % Express g0(x)+g(x)*d=0 as Ad*d = b
    if any(isdec)
        [z,b,Aq,Ad] = gramconstraint(g0,g);
        Ad = -Ad;
    else
        [z,b,Aq] = gramconstraint(g0);
    end
    Nconstr = length(b);
    
    % Stack A*d=b into As*xs=bs where the free dec. vars fv precede
    % the slack variables, y, which precede the SOS dec vars sv.
    tmp = sparse(Nconstr,Nfv+Nlpv+Nsosv);
    tmp(:,dv2xs(didx)) = Ad;
    As = [As; tmp];
    bs = [bs; b];
end

% Define LP constraints: A*d+y=b, y>=0
K.l = Nclp;
if K.l>0
    s = standform.linprog.data;
    A = s{1};
    b = s{2};
    d = s{3};
    
    % Check that constraint is linear in the decision variables
    if ispvar(d)
        % Find the decision variables involved in the LP constraint
        [isdec,didx]=ismember(char(d),decvars);
    elseif isempty(d)
        isdec=0;
        didx=[];
    else
        % Error: Nonlinear dependence on decision variables
        tmpidx = find( sum(d.degmat,2) > 1);
        dterm = polynomial(1,d.degmat(tmpidx(1),:),d.varname,[1 1]);
        
        didx = find( d.coef(tmpidx(1),:) );
        ridx = find( A(:,didx(1)) );
        error( LOCALerrstr(dterm,standform.linprog.idx(ridx)) )
    end
    
    % Stack A*d+y=b into As*xs=bs where the free dec. vars fv precede
    % the slack variables, y, which precede the SOS dec vars sv.
    tmp = sparse(Nclp,Nfv+Nlpv+Nsosv);
    tmp(:,dv2xs(didx)) = A;
    tmp(:,Nfv+(1:Nclp)) = speye(Nclp);
    As = [As; tmp];
    bs = [bs; b];
end

% Define SOS Variable Constraints: M>=0
zs = cell(Ncsosdv,1);
K.s = [];
for i1=1:Ncsosdv
    s = standform.sosdecvar.data(i1,:);    
    zs{i1} = s{1};
    K.s = [K.s length(s{1})];
end

% Define SOS Inequality Constraints: Ad*d + Aq*Q(:) = b and Q>=0
zi = cell(Ncsos,1);
ptr = Nfv+Nlpv+Nsosv;
for i1=1:Ncsos
    s = standform.sos.data(i1);
    
    % Collect s(x,d) into the form g0(x)+g1(x)*d1 + ... + gN(x)*dN
    % where di are decision variables.
    [g0,g,d] = collect(s,x);
    if ispvar(d)
        % Find the decision variables involved in the SOS constraint
        [isdec,didx]=ismember(char(d),decvars);
    elseif isempty(d)
        isdec=0;
        didx=[];
    else
        tmpidx = find( sum(d.degmat,2) > 1);
        dterm = polynomial(1,d.degmat(tmpidx(1),:),d.varname,[1 1]);
        error( LOCALerrstr(dterm,standform.sos.idx(i1)) )
    end
    
    % Express g0(x)+g(x)*d is SOS as
    %     g0(x)+g(x)*d = z(x)'*Q*z(x),  Q>=0
    % Then convert the polynomial equality constraint into:
    %           Ad*d + Aq*Q(:) = b
    if any(isdec)
        [z,b,Aq,Ad] = gramconstraint(g0,g);
        Ad = -Ad;
    else
        [z,b,Aq] = gramconstraint(g0);
    end
    zi{i1} = z;
    Nconstr = length(b);
    
    % Stack Ad*d + Aq*Q(:) = b into As*xs=bs where the psd slack vars,
    % Q(:), follow the free dec vars fv, the LP slack vars y, and the
    % SOS dec vars sv.
    K.s = [K.s length(z)];
    tmp = [sparse(Nconstr,ptr) Aq];
    if any(isdec)
        tmp(:,dv2xs(didx)) = Ad;
    end
    if isempty(As)
        As = tmp;
    else
        As = [ [As sparse(size(As,1),size(Aq,2))] ; tmp];
    end
    bs = [bs; b];
    ptr = size(As,2);
end
z = [zs;zi];

%---------------------------------------------------------------------
% Monomial and decision variable reduction
%---------------------------------------------------------------------
if strcmpi(opts.simplify,'on');
    [As,bs,K,z,dv2xs,Nfv,feas,zrem] = ...
        sospsimplify(As,bs,K,z,dv2xs,Ncsosdv);
else
    zrem = cell(length(K.s),1);
    feas = 1;
end

%---------------------------------------------------------------------
% Create the objective function vector
%---------------------------------------------------------------------
if isempty(obj)
    % Feasibility Problem
    cs = sparse( size(As,2) , 1 );
else    
    % Objective function is obj = cs'*xs
    [g0,g,d] = collect(obj,x);
    if (g0.nvar ~=0) || (g.nvar~=0)
        % Objective depends on polynomial vars
        gtmp = [g0(:);g(:)];
        tmpidx = find( sum(gtmp.degmat,2) > 0);
        gterm = polynomial(1,gtmp.degmat(tmpidx(1),:),gtmp.varname,[1 1]);
        error( LOCALerrstr(gterm,-2) )
    elseif ispvar(d)
        % Find the decision variables involved in the objective function
        [isdec,didx]=ismember(char(d),decvars);
    elseif isempty(d)
        %isdec=0;
        didx=[];
    else
        % Objective has nonlinear dependence on dec vars
        tmpidx = find( sum(d.degmat,2) > 1);
        dterm = polynomial(1,d.degmat(tmpidx(1),:),d.varname,[1 1]);
        error( LOCALerrstr(dterm,-1) )
    end
    %g0 = double(g0);
    
    % Create temporary indices for removed variables
    Nv = size(As,2);
    idx = find(dv2xs==0);
    Nrem = length(idx);
    tmpdv2xs = dv2xs;
    tmpdv2xs(idx) = Nv+(1:Nrem);
    
    % Create obj function and strip away indices for removed vars
    cs = sparse( Nv+Nrem , 1 );
    cs(tmpdv2xs(didx)) = double(g);
    cs = cs(1:Nv);
end

%---------------------------------------------------------------------
% Create info data structure
%---------------------------------------------------------------------

% sosdata
sosdata.decvars = decvars;

% sdpdata
sdpdata.A = As;
sdpdata.b = bs;
sdpdata.c = cs;
sdpdata.K = K;

% sos2sdp: z, dv2x, dv2y, Nfv, Nlpv, orderidx
sos2sdp.z = z;
sos2sdp.zremoved = zrem;
sos2sdp.dv2x = dv2xs;
sos2sdp.dv2y = [];

numvar.free = Nfv;
numvar.lpslack = Nlpv;
numvar.sosvar = Nsosv;
numvar.sosgram = size(As,2)-(Nfv+Nlpv+Nsosv);

sos2sdp.numvar = numvar;

% Pack all data into info
info.sosdata = sosdata;
info.sdpdata = sdpdata;
info.sos2sdp = sos2sdp;
info.opts = opts;

%---------------------------------------------------------------------
% Solve the problem
%---------------------------------------------------------------------
if strcmpi(opts.solver,'setup')
    % Return with problem formulation for setup call
    sdpsol.x = [];
    sdpsol.y = [];
    sdpsol.solverinfo = [];
    info.sdpsol = sdpsol;
    info.obj = [];
    
    dopt =[];
    sossol = [];
    return;
end
if feas
    if ~isempty(As)
        [xs,ys,solverinfo] = solvesossdp(sdpdata,opts);
    else
        % 6/22/10: Abhijit found an example that was feasible but
        % SOSSIMPLIFY removed all variables.
        % XXX This is only the correct solverinfo for Sedumi
        xs = zeros(size(As,2),1);
        ys = zeros(size(As,1),1);
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
    xs = zeros(size(As,2),1);
    ys = zeros(size(As,1),1);
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
info.sdpsol = sdpsol;

%---------------------------------------------------------------------
% Create sossol and dopt for output
%---------------------------------------------------------------------

% Create dopt
if Ndv ~= 0
    dvals = zeros(Ndv,1);
    idx = find(dv2xs~=0);
    dvals( idx ) = xs( dv2xs(idx) );
    %dopt = [decvars(:) num2cell(dvals)];
    dopt = [polynomial(decvars(:)) dvals];
else
    dopt = [];
end

% Fill in the minimal objective value
if isempty(obj)
    info.obj = [];
elseif Ndv~=0
    info.obj = double(subs(obj,dopt));
else
    info.obj = double(obj);
end

% Create sossol
sossol = [];
if nargout==3 
    cnt = Nfv;
    
    for i1=1:Np
        % Evaluate pconstr at dopt
        p = pconstr(i1);
        if Ndv~=0
            p = subs(p,dopt);
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

% LOCAL function to generate error string for nonlinear dependence
% on decision variables
function estr = LOCALerrstr(dterm,idx)

if idx==-1
    % Objective function has nonlinear dependence on dec vars
    dstr = char(dterm);    
    errstr1 = 'Objective function contains a ';
    errstr2 = ['term that depends on ' dstr{1} '.'];
    errstr3 = 'Objective function must be linear in the decision variables.';
    estr = sprintf([errstr1 errstr2 '\n' errstr3]);
elseif idx==-2
    % Objective function depends on poly vars
    dstr = char(dterm);    
    errstr1 = 'Objective function contains a ';
    errstr2 = ['term that depends on ' dstr{1} '.'];
    errstr3 = 'Objective function can not depend on polynomial independent variables.';
    estr = sprintf([errstr1 errstr2 '\n' errstr3]);
else
    % Nonlinear dependence in constraint
    dstr = char(dterm);
    
    errstr1 = ['Constraint ' int2str(idx) ' contains a '];
    errstr2 = ['term that depends on ' dstr{1} '.'];
    errstr3 = 'Constraints must be linear in the decision variables.';
    estr = sprintf([errstr1 errstr2 '\n' errstr3]);
end