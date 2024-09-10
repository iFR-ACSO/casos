function [vol,volstd] = pvolume(casospoly,v,domain,npts)
% DESCRIPTION 
%   Estimate the volume contained in the set {x : p(x)<=v} using Monte
%   Carlo sampling.  npts are drawn uniformly from a hypercube and the
%   number of points, nin, contained in the set { x : p(x) <= v} is
%   counted.  The volume is estimated as vol = nin/npts.  An estimate
%   of the standard deviation of this volume is also computed.
%   See Ref: Monte Carlo Integration Wikipedia Entry
%
% INPUTS 
%   p: 1-by-1 polynomial of n variables (casospoly)
%   v: scalar specifying the sublevel of the polynomial (Default: v=1)
%   domain: n-by-3 array specifying the sampling hybercube.  domain(i,1)
%           is a indeterminate in p and domain(i,2:3) specifies the min and max
%           values of the cube along the specified variable direction,
%              [X1, X1min, X1max; ...; Xn, Xnmin, Xnmax]
%         (Default: domain = [-1 1] along all variable directions)
%   npts: scalar specifying the number of sample points 
%          (Default: npts = 1e4)
%
% OUTPUTS  
%   vol: Volume estimate of { x : p(x)<= v}
%   volstd: Standard deviation of the volume estimate.
%  
% SYNTAX 
%  pvolume(p)
%  pvolume(p,v)
%  pvolume(p,v,domain)
%  pvolume(p,v,domain,npts)
%  [vol,volstd] = pvolume(p,v,domain,npts)
%
% EXAMPLE
%   x = casos.Indeterminates('x', 3, 1);
%   V = 1*(x'*x);
%   npts = 1e5;
%   domain = [x(1) -2 2; x(2) -1 1; x(3) -1 1]; 
%   [vol,stdvol] = pvolume(V, 1, domain, npts)

if nargin==1
    v=[];
    domain=[];
    npts=[];
elseif nargin==2    
    domain=[];
    npts=[];
elseif nargin==3
    npts=[];
end

% Default contour
if isempty(v)
    v=1;
end

% Default npts
if isempty(npts)    
    npts = 1e4;
end


if isempty(domain)
    nvar = casospoly.nvars;
    xmin = -ones(1,nvar);
    xmax = ones(1,nvar);
else
    nvar = size(domain,1);
    
    L.type = '()';
    L.subs = {':',1};
    var = subsref(domain,L);
    L.subs = {':',2};
    xmin = full( subsref(domain,L) )';
    L.subs = {':',3};
    xmax = full( subsref(domain,L) )';

    % obtain permutation matrix
    permutation = full(nabla(1.*casospoly.indeterminates, var));
    if rank(permutation)~= nvar
        error('Wrong domain!')
    end
    % reorder
    xmin = xmin*permutation';
    xmax = xmax*permutation';
end    

% Hypercube volume 
xdiff = xmax-xmin;
cubevol = prod(xdiff);

% Generate samples
xdiff = repmat(xdiff,[npts 1]);
xmin = repmat(xmin,[npts 1]);
xpts = xmin+xdiff.*rand(npts,nvar);

% Convert matrix columns to a cell array, each cell contains one column
columnCells = mat2cell(xpts, size(xpts, 1), ones(1, size(xpts, 2)));

% from casos poly to casadi function
polyFun = to_function(casospoly);

% evaluate casadi function
pv = full(polyFun(columnCells{:}));


% Estimate Volume and Standard Deviation
% Ref: Monte Carlo Integration Wikipedia Entry
f = pv<=v;
fmean = sum(f)/npts;
fvar = sum( (f-fmean).^2 )/(npts-1);

vol = fmean*cubevol;
volvar = cubevol^2*fvar/npts;
volstd = sqrt(volvar);

