function [gbnds,sout,sossol,info]=pcontain(p1,p2,z,opts)
% function [gbnds,s,sossol]=pcontain(p1,p2,z,opts)
%
% DESCRIPTION 
%   This function maximizes g subject to the set containment constraint:
%           { x : p2(x)<= g } is a subset of { x : p1(x)<= 0 }
%   The set containment constraint is relaxed using a generalization of
%   the S-procedure, i.e. the set containment holds if there exists a 
%   function s(x) such that s and -( p1+(g-p2)*s ) are sums of squares.
%   This solved as the following SOS problem:
%        max g
%        subject to: s in SOS,  -( p1+(g-p2)*s ) in SOS
%   Bisection is required to find the optimal g since the constraint is,
%   in general, bilinear in g and s(x). The bisection stops when the 
%   upper/lower bounds on g are within a specified absolute tolerance.
%
% INPUTS 
%   p1,p2: Polynomials used to describe the sets in the containment.
%   z: Nz-by-1 column vector of monomials used to specify the 
%      SOS decision variable s(x) in the Gram matrix form, 
%      s(x)=z(x)'*C*z(x). [Optional] Default is z=monomials(x,[mind,maxd])
%      where maxd is chosen such that maxdeg(p2*s) >= maxdeg(p1) 
%      and mind is chosen such that s.mindeg>=p1.mindeg.
%   opts: Options for optimization.  See GSOSOPTIONS for more details.
%     The minobj and maxobj options refer to known bounds on the
%     maximum value of gamma.
%     
% OUTPUTS
%   gbnds: 1-by-2 vector [glb, gub] providing a lower bound glb and
%      and upper bound gub on the maximum value of g. gbnds will be
%      empty if no feasible lower bound is found.
%   s: Multplier function proving set containment with glb.  s is SOS
%      and satisfies -(p1+(glb-p2)*s) in SOS. s will be empty if no
%      feasible lower bound is found.
%   sossol: 2-by-1 structure array with fields p, z, and Q.
%      sossol(1).p is the optimal multiplier s and sossol(2).p is the set 
%      containment polynomial -(p1+(glb-p2)*s ) at the optimal s.  
%      sossol(i).z is a vector of monomials in x and sossol(i).Q is a
%      positive semidefinite matrix such that p=z'*Q*z (i=1,2). This is a 
%      Gram matrix decomposition which proves the two constraints are SOS.
%   info: Information structure returned by gsosopt.
%
% SYNTAX
%   [gbnds,s,sossol,info]=pcontain(p1,p2)
%   [gbnds,s,sossol,info]=pcontain(p1,p2,z)
%   [gbnds,s,sossol,info]=pcontain(p1,p2,z,opts)
%
% EXAMPLE
%   % Maximize size of an ellipsoid contained in unit disk
%   pvar x1 x2;
%   x = [x1;x2];
%   p1 = x'*x-1;
%   p2 = 4*x1^2+9*x2^2;
%   z = 1;
%   [gbnds,sout,sossol,info]=pcontain(p1,p2,z)
%
% See also gsosopt, gsosoptions

%  4/27/2009 PJS  Initial coding based on getgamma and getbeta
%  12/7/2009 PJS  Update to call gsosopt
%  11/01/10  PJS  Update for polyconstr, sosoptions, == constraints


% XXX-- (?) Add alternative formulation for an outer set approx:
%   min beta s.t. {p2<=0} contained in {p1<=beta}

% XXX-- Add exact solution for the case where p1, p2 are quadratic
% SOS forms?  This can be solved by computing a max gen. eigenval.
% 
%  max g s.t. {x'*Q*x<=g} <= {x'*P*x<=1}
% 
% Solution is:
%  a) lmax = max(eig(P,Q))
%  b) gmax = 1/lmax
% 
% % Define Random Ellipsoids
% X = mpvar('x',2,1);
% P = randn(2); P = P*P'+0.1*eye(2);
% Q = randn(2); Q = Q*Q'+0.1*eye(2);
% 
% % Solve generalized eigenvalue problem
% [EVEC,EVAL] = eig(P,Q);
% [lmax,idx] = max(diag(EVAL));
% v = EVEC(:,idx);
% gmax = 1/lmax;
% 
% % Plot contours 
% dom = [-2 2 -2 2];
% pcontour( X'*Q*X , gmax, dom, 'b'); hold on;
% pcontour( X'*P*X , 1, dom, 'r'); 
% 
% % Plot location where ellipsoids touch
% c = v'*P*v;
% v = v/sqrt(c);
% plot(v(1),v(2),'ko'); hold off;

% XXX-- Document info returned from gsosopt?

%------------------------------------------------------------------
% Error Checking
%------------------------------------------------------------------

%------------------------------------------------------------------
% Get options or set defaults for unspecified options
%------------------------------------------------------------------
if nargin==2
    z = [];
    opts = [];
elseif nargin==3
    opts = [];
end

%------------------------------------------------------------------
% Initialize variables
%------------------------------------------------------------------

% Get polynomial variables
x = unique([p1.varname; p2.varname]);

% Set default options
if isempty(opts)
    opts = gsosoptions;
end

% Set default multiplier
if isempty(z)
    % This is roughly based on guidelines from Appendix A.1.1 
    % of the Ph.D. Thesis by Weehong Tan
    zmind = ceil(p1.mindeg/2);
    zmaxd = ceil( (p1.maxdeg-p2.maxdeg) / 2 );    
    if zmaxd<zmind
        z=1;
    else
        z = monomials(x, zmind:zmaxd );
    end
end

% Call GSOSOPT to solve:
%     max g
%     s in SOS 
%     -( p1+(g-p2)*s ) = -(p1+g*s-p2*s) in SOS
% GSOSOPT solves a minimization using t:=-g
%     min t
%     s in SOS
%     (p2*s-p1) + t*s in SOS
t = pvar('g');
s = sosdecvar('c',z);

sosc(1) = s>=0;
sosc(2) = (p2*s-p1)+t*s>=0;

gopts = opts;
gopts.minobj = -opts.maxobj;
gopts.maxobj = -opts.minobj;
if strcmp(opts.display,'on')
    gopts.display = 'pcontain'; % Undocumented syntax for proper display
end

if nargout>=3
    [info,dopt,sossol] = gsosopt(sosc,x,t,gopts);
else
    [info,dopt] = gsosopt(sosc,x,t,gopts);
end

if ~isempty(info.tbnds)
    gbnds = [-info.tbnds(2) -info.tbnds(1)];    
    sout = subs(s,dopt);
else
    gbnds = [];
    sout = polynomial([]);
end

