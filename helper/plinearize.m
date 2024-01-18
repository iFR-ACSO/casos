function [A,B,f0] = plinearize(f,x,u,x0,u0)



% Error Checking
if nargin==2
    %   [A,f0] = plinearize(f,x)
    u = [];
    x0 = [];
    u0 = [];
elseif nargin==3

        x0 = [];
        u0 = [];        

elseif nargin==4
    %   [A,B,f0] = plinearize(f,x,u,x0)
    u0 = [];
end

% Default trim values
if isempty(x0)
    Nx = length(x);
    x0 = zeros(Nx,1);
end
if isempty(u0)
    Nu = length(u);
    u0 = zeros(Nu,1);
end
   
% Evaluate function at trim
v = [x(:); u(:)];
v0 = [x0(:); u0(:)];

f0 = subs(f,v,v0);
f0 = double(f0);

% State matrix 
dfdx = nabla(f,x);
A = subs(dfdx,v,v0);
A = double(A);

% Input matrix
if ~isempty(u)
    dfdu = nabla(f,u);
    B = subs(dfdu,v,v0);
    B = double(B);
else
    B=f0;
end
