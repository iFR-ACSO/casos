function varargout = pcontour3(poly,lvlset,domain,noPoints,varName)
% Plot 2D contour plot of a polynomial.
%
% Input:
%   poly       - casos polynomial \in \mathbb R[x1,x2,x3]
%   lvlset     - desired level set    (double/integer)                  
%                default: 1
%   domain     - domain for plotting                                    
%                default: [-1 1 -1 1 -1 1]
%   noPoints   - number of discrete points for level set interpolation  
%                default: [100 100]
%   var        - varible names (char)
%
% Output: 
%   none or
%   C - Contour matrix 
%   h - contour object (can be used to set properties of contour plot)
%
% Example:  Plot zero level set of an ellipse on the domain 
%           [-2 2] \times [-2 2] and color line as black dashed.
%
%   x = casos.PS('x',3,1);
%   poly = x(1)^2 + x(2)^2 + x(3)^2 - 1;
%        
%   [C,h] = pcontour3(poly,0,[-2 2 -2 2 -2 2],'k--');
%
%   pcontour(poly,0,[-2 2 -2 2],'k--');
%


%% set default values (similar to sosopt)
if nargin==1
    lvlset=[];
    domain=[];
    noPoints=[];
    varName=[];
elseif nargin==2    
    domain=[];
    noPoints=[];
    varName=[];
elseif nargin==3
    noPoints=[];
    varName=[];
elseif nargin==4
    noPoints=[];
    varName=[];
elseif nargin==5
    varName=[];
end

% Default contour
if isempty(lvlset)
    lvlset=1;
end


% Default noPoints
if isempty(noPoints)    
    Nx = 100;
    Ny = 100;
    Nz = 100;
else
    Nx = noPoints(1);
    if length(noPoints)==1
        Ny = Nx;
    else
        Ny = noPoints(2);
    end
end

% Get variable names
if isempty(varName)
    x = poly.indeterminates.str{1};
    y = poly.indeterminates.str{2};
    z = poly.indeterminates.str{3};
elseif ispvar(varName) || iscellstr(varName)
    if ispvar(varName)
        varName = char(varName);
    end
    x = varName{1};
    y = varName{2};
    z = varName{3};
else
    error('varName must be a vector of pvars');
end


if isempty(domain)
    domain = [-1 1 -1 1 -1 1];
end    
    
%% Plot contour

% generate grid
xgrid = linspace(domain(1),domain(2),Nx);
ygrid = linspace(domain(3),domain(4),Ny);
zgrid = linspace(domain(5),domain(6),Nz);
[xgrid,ygrid,zgrid] = meshgrid(xgrid,ygrid,zgrid);


% from casos poly to casadi function
polyFun = to_function(poly);

% evaluate casadi function
gridval = full(polyFun([xgrid(:)';ygrid(:)';zgrid(:)']));

% reshape to grid
gridval = reshape(gridval,size(xgrid));


% execute matlab isosurface function
if nargout==0
    isosurface(xgrid,ygrid,zgrid,gridval,lvlset);    
    axis(domain)
    xlabel(x)
    ylabel(y)
    zlabel(z)
elseif nargout==1
    F=isosurface(xgrid,ygrid,zgrid,gridval,lvlset);    
    varargout = {F};
elseif nargout==2
    [F,V]=isosurface(xgrid,ygrid,zgrid,gridval,lvlset);    
    varargout = {F,V};
elseif nargout==3
    [F,V,C]=isosurface(xgrid,ygrid,zgrid,gridval,lvlset);    
    varargout = {F,V,C};
end


end