function [vol,volstd] = pvolume(p,v,domain,npts)
% DESCRIPTION 
%   Estimate the volume contained in the set {x : p(x)<=v} using Monte
%   Carlo sampling.  npts are drawn uniformly from a hypercube and the
%   number of points, nin, contained in the set { x : p(x) <= v} is
%   counted.  The volume is estimated as vol = nin/npts.  An estimate
%   of the standard deviation of this volume is also computed.
%   See Ref: Monte Carlo Integration Wikipedia Entry
%
% INPUTS 
%   p: 1-by-1 polynomial of n variables 
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
    v = 1;
elseif ~isa(v, 'double') 
    assert(v.is_constant, 'Level set value must be a constant');
    v = full(casos.PD(v));
end

% Default npts
if isempty(npts)    
    npts = 1e4;
end

if isempty(domain)
    % Default domain: [-1, 1] for each variable
    nvar = p.nvars;
    xmin = -ones(1,nvar);
    xmax = ones(1,nvar);
else
    % Extract variables and bounds from the domain
    nvar = size(domain,1); 
    var  = domain(:,1);                             % First column: variables
    bounds = horzcat(domain(:,2), domain(:,3));     % Get domain bounds

    % Check whether an object can be converted to casos.Indeterminates
    assert(is_indet(var), 'First column of the domain must be a vector of indeterminate variables.');

    % Check if the domain does not contain symbolic variables
    assert(bounds.is_constant, 'Domain must not contain symbolic variables.');

    % Convert bounds to PD followed by conversion to double via full
    bounds = full(casos.PD(bounds));

    [xmin, xmax] = deal(bounds(:,1), bounds(:,2));

    % Convert variables to indeterminates (preserves sorting)
    indets = casos.Indeterminates(var);

    % Match polynomial indeterminates with domain variables
    [tf, loc] = ismember(p.indeterminates, indets);
    assert(all(tf), 'Invalid domain.')

    % Reorder bounds to match polynomial indeterminates
    xmin = xmin(loc);
    xmax = xmax(loc);
end    

% Hypercube volume 
xdiff = xmax-xmin;
cubevol = prod(xdiff);

% Generate samples
xdiff = repmat(xdiff,[1 npts]);
xmin = repmat(xmin,[1 npts]);
xpts = (xmin+xdiff.*rand(nvar, npts))';

% Convert matrix columns to a cell array, each cell contains one column
columnCells = mat2cell(xpts, size(xpts, 1), ones(1, size(xpts, 2)));

% From casos poly to casadi function
polyFun = to_function(p);

% Evaluate casadi function
pv = full(polyFun(columnCells{:}));

% Estimate Volume and Standard Deviation
% Ref: Monte Carlo Integration Wikipedia Entry
f = pv<=v;
fmean = sum(f)/npts;
fvar = sum( (f-fmean).^2 )/(npts-1);

vol = fmean*cubevol;
volvar = cubevol^2*fvar/npts;
volstd = sqrt(volvar);

