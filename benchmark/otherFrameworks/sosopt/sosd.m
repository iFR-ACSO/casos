function  [info,dopt,sossol] = sosd(sosdata,opts)
% function [info,dopt,sossol] = sosd(sosdata,opts)
%
% DESCRIPTION
%   This function solves the SOS optimization or feasibility problem in
%   kernel (dual) form.  An SOS optimization contains constraints of the
%   form (see sosopt for more details):
%               s(x,d) is sum of squares
%   where d are the decision variables and x are the polynomial variables.
%   Each SOS constraint is equivalent to:
%               s(x,d)=z(x)'*Q*z(x)
%               Q>=0
%   The constraint s(x,d)=z(x)'*Q*z(x) is simply a linear equality
%   constraint between the entries of Q, the decision variables d, and
%   the numeric coefficients of s.
%
%   This function solves the optimization by constructing a particular
%   solution and a basis for the null space of this linear equality
%   constraint.  Specifically, matrices Q0, D1,...,Dk are constructed
%   such that z(x)'*(Q0 - sum_i d(i)*Di)*z(x) is equal to s(x,d).  Also,
%   a set of basis matrices N1,...,Nj are constructed such that
%   z(x)'*Ni*z(x) = 0.  Finally, the SOS constraint is replaced with
%   the linear matrix inequality constraint:
%            Q0 - sum_i di*Di -sum_i ni *Ni >=0.
%
% See also sosopt, sosoptions, sosp

% 2/29/08  PJS     Initial Coding
% 11/11/08 PJS     Simple mod to handle sos var representation
% 4/25/09 PJS      Clean up documentation and code for pure sos vars
% 5/6/09 PJS       Modify to allow set-up call
% 11/22/10 PJS     Update for polyconstr, sosoptions, == constraints

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
Ndv = length(decvars);

if Nceq>0
    % XXX TODO Incorporate equality constraints
    error(['Kernel form currently does not handle equality '...
        'constraints. Use image form instead.']);
end

%---------------------------------------------------------------------
% Define SOS/LP Cost and Constraints in Sedumi Dual Form
%       max bs'*ys   s.t. cs-As'*ys  >= 0
%
% SOS Constraints: The constraint s(x,d) is SOS is equivalent to the
%     constraints s(x,d) = z(x)'*Q*z(x) and Q>=0.  These can be expressed
%     as Q0 - sum_i di*Di -sum_i ni *Ni >=0.
% LP Constraints: A*d<=b is defined as b-A*d>=0.
%
% Sedumi dual form orders the LP constraints before the LMI constraints.
% The Sedumi decision variables are ordered as:
%       xs=[d; n1; ....; nlast]
% where d are the free decision variables, ni are the null space
% variables needed to express the SOS constraints.
%---------------------------------------------------------------------

% Define LP constraints: b-A*d>=0
K.l = Nclp;
Ast = sparse([]);    % Ast = As'
cs = sparse([]);
if K.l>0
    % Get LP constraint data
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
    
    % Stack LP constraints below preceding constraints
    tmp = sparse(Nclp,Ndv);
    tmp(:,didx) = A;
    if isempty(Ast)
        Ast = tmp;
    else
        Ast = [Ast; tmp];
    end
    cs = [cs; b];
end

% Define SOS Constraints: Q0 - sum_i di*Di -sum_i ni *Ni >=0
K.s = [];
zi = cell(Ncsos,1);
for i1=1:(Ncsosdv+Ncsos)
    % i^th constraint is s(x,d) is SOS where d are the decision vars
    if i1<=Ncsosdv
        % Pure SOS Dec Var constraint
        data = standform.sosdecvar.data(i1,:);
        z = data{1};
        s = data{3};
    else
        % Generic SOS Constraint
        s = standform.sos.data(i1-Ncsosdv);
        z = [];
    end
    
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
    %     g0(x)+g(x)*d = z(x)'*Q*z(x),  Q>=0  (*)
    % Then find Q0, Ni and Di such that:
    %     z'*Q0*z = g0, z'*Ni*z=0, and z'*Di*z = gi
    % Thus (*) is equivalent to the LMI constraint:
    %     Q0 + sum_i ni*Ni + sum_i di*Di>=0
    if any(isdec)
        [z,Q0,N,D] = gramsol(g0,g,z);
    else
        [z,Q0,N] = gramsol(g0,[],z);
    end
    Q0 = Q0(:);
    zi{i1} = z;
    
    % Add this constraint to the SOS problem
    % Null space variables are added to the end of the Sedumi decision
    % variable vector.  The minus signs account for the sign in the
    % Sedumi dual constraint: c-A'*y>=0
    lz = length(z);
    K.s = [K.s lz];
    if all(size(Ast)==[0 0]) %isempty(Ast)
        Ast = [ sparse(lz^2,Ndv) -N];
    else
        Ast = blkdiag(Ast,-N);
    end
    if any(isdec)
        Ast(end-lz^2+1:end,didx) = -D;
    end
    cs = [cs; Q0];
end
As = Ast';
z = zi;

%---------------------------------------------------------------------
% Create the objective function vector
%---------------------------------------------------------------------
if isempty(obj)
    % Feasibility Problem
    bs = sparse( size(Ast,2) , 1 );
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
    
    % SOS programming is formulated as a minimization but the Sedumi
    % dual form is a maximization.  Flip sign of bs to account for this.
    bs = sparse( size(Ast,2) , 1 );
    bs(didx) = -double(g);
end

%---------------------------------------------------------------------
% Create info data structure
%---------------------------------------------------------------------

% sdpdata
sdpdata.A = As;
sdpdata.b = bs;
sdpdata.c = cs;
sdpdata.K = K;

% sos2sdp: z, dv2x, dv2y, Nfv, Nlv, orderidx
sos2sdp.z = z;
sos2sdp.dv2x = [];
sos2sdp.dv2y = 1:Ndv;
sos2sdp.Nfv = [];
sos2sdp.Nlv = [];
sos2sdp.orderidx = [];

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
if isempty(bs)
    % 11/29/10 -- sosopt(x^2) in kernel form gives an SDP with
    % both A and b as 0-by-1. This causes Sedumi to error out.
    % XXX This is only the correct solverinfo for Sedumi
    e = eigK(full(sdpdata.c),K);
    tol = 1e-8;  % XXX
    xs = zeros(length(sdpdata.c));
    ys = zeros(0,1);
    if min(real(e))>= -tol
        solverinfo.feasratio = 1;
        solverinfo.pinf = 0;
        solverinfo.dinf = 0;
        solverinfo.numerr = 0;
    else
        solverinfo.feasratio = 1;
        solverinfo.pinf = 1;
        solverinfo.dinf = 1;
        solverinfo.numerr = 0;
    end
        
else
    [xs,ys,solverinfo] = solvesossdp(sdpdata,opts);
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
    %dopt = [decvars(:) num2cell(ys(1:Ndv))];
    dopt = [polynomial(decvars(:)) ys(1:Ndv)];
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
    for i1=1:Np
        % Evaluate pconstr at dopt
        p = pconstr(i1);
        if Ndv~=0
            p = subs(p,dopt);
        end
        sossol(i1,1).p = p;
    end
    
    %     % Fill sossol for equality constraints
    %     % XXX Return anything for equality constraints?
    %     idx=standform.equality.idx;
    %     for i1=1:Nceq
    %         sossol(idx(i1)).z = [];
    %         sossol(idx(i1)).Q = [];
    %     end
    
    % Fill sossol for LP constraints
    idx = standform.linprog.idx;
    if Nclp>0
        Qvec = cs(1:Nclp) - Ast(1:Nclp,:)*ys;
        for i1=1:Nclp
            % Unscale: XXX Currently no LP scaling.
            sossol(idx(i1)).z = 1;
            sossol(idx(i1)).Q = Qvec(i1);
        end
    end
    
    % The call to subs is only moderately fast so only do this if needed
    cnt = K.l;
    for i1=1:(Ncsosdv+Ncsos)
        if i1<=Ncsosdv
            % Pure SOS Dec Var constraint
            idx=standform.sosdecvar.idx(i1);
            scl = 1;
        else
            % Generic SOS Constraint
            idx=standform.sos.idx(i1-Ncsosdv);
            scl = scaling.sos(i1-Ncsosdv);
        end
        xidx = cnt+(1:K.s(i1)^2);
        Qvec = cs(xidx) - Ast(xidx,:)*ys;
        Q = reshape( Qvec , K.s(i1), K.s(i1) );
        
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

