function varargout = pcontour(poly,lvlset,domain,linespec,noPoints,varName)
% Plot 2D contour plot of a polynomial.
%
% Input:
%   poly       - casos polynomial \in \mathbb R[x1,x2]
%   lvlset     - desired level set    (double/integer)                  
%                default: 1
%   domain     - domain for plotting                                    
%                default: [-1 1 -1 1]
%   linespec   - defined line color and type for level set              
%                default: 'b'
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
%   x = casos.PS('x',2,1);
%   poly = x(1)^2 + x(2)^2 - 1;
%        
%   [C,h] = pcontour(poly,0,[-2 2 -2 2],'k--');
%
%   pcontour(poly,0,[-2 2 -2 2],'k--');
%


%% set default values (similar to sosopt)
if nargin==1
    lvlset=[];
    domain=[];
    linespec=[];
    noPoints=[];
    varName=[];
elseif nargin==2    
    domain=[];
    linespec=[];
    noPoints=[];
    varName=[];
elseif nargin==3
    linespec=[];
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


% Default linespec
if isempty(linespec)
    linespec = 'b';
end

% Default noPoints
if isempty(noPoints)    
    Nx = 100;
    Ny = 100;
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
elseif ispvar(varName) || iscellstr(varName)
    if ispvar(varName)
        varName = char(varName);
    end
    x = varName{1};
    y = varName{2};
else
    error('varName must be a vector of pvars');
end


if isempty(domain)
    domain = [-1 1 -1 1];
end    
    
%% Plot contour

% generate grid
xg = linspace(domain(1),domain(2),Nx);
yg = linspace(domain(3),domain(4),Ny);
[xg,yg] = meshgrid(xg,yg);


% from casos poly to casadi function
polyFun = to_function(poly);

% evaluate casadi function
pgrid = full(polyFun(xg,yg));

% reshape to grid
pgrid = reshape(pgrid,size(xg));

if length(lvlset) == 1
    % Single contour syntax for contour function
    lvlset = [lvlset lvlset];
end

% execute matlab contour plot function
if nargout==0
    contour(xg,yg,pgrid,lvlset,linespec);    
else
    [C,h]=contour(xg,yg,pgrid,lvlset,linespec);     
    varargout = {C,h};
end
xlabel(x)
ylabel(y)


end