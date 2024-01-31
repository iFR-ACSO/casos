function [vol,volstd] = pvolume(casospoly,v,domain,npts)

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
    var  = casospoly.indeterminates.str; 
    nvar = length(var);
    xmin = -ones(1,nvar);
    xmax = ones(1,nvar);
else
    nvar = size(domain,1);
    
    L.type = '()';
    L.subs = {':',1};
    var = subsref(domain,L);
    L.subs = {':',2};
    xmin = double( subsref(domain,L) )';
    L.subs = {':',3};
    xmax = double( subsref(domain,L) )';
end    

% Hypercube volume 
xdiff = xmax-xmin;
cubevol = prod(xdiff);

% Generate samples
xdiff = repmat(xdiff,[npts 1]);
xmin = repmat(xmin,[npts 1]);
xpts = xmin+xdiff.*rand(npts,nvar);

% from casos poly to casadi function
polyFun = to_function(casospoly);

% evaluate casadi function
pv = full(polyFun(xpts'));


% Estimate Volume and Standard Deviation
% Ref: Monte Carlo Integration Wikipedia Entry
f = pv<=v;
fmean = sum(f)/npts;
fvar = sum( (f-fmean).^2 )/(npts-1);

vol = fmean*cubevol;
volvar = cubevol^2*fvar/npts;
volstd = sqrt(volvar);

