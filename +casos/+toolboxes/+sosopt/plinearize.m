function [A,B,f0] = plinearize(f,x,u,x0,u0)
% Linearize polynomial function.

Nx = length(x);
Nu = length(u);

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
    x0 = zeros(Nx,1);
end
if isempty(u0)
    u0 = zeros(Nu,1);
end
   
% Evaluate function at trim
f0 = full(subs(f,[x;u],[x0(:); u0(:)]));

% State matrix 
fx = nabla(f,x);

% check if any state is member of A matrix; if yes evaluate it at trim point; 
% otherwise convert to double
% if any(ismember(fx.indeterminates.str,x.str))

    A = casos.PD(subs(fx,[x;u],[x0;u0]));

% else
%     A = full(fx);
% end

% Control matrix
if ~isempty(u)
    fu = nabla(f,u);
    
    % check if any control is member of B matrix; if yes 
    % evaluate it at trim point; otherwise convert to double
    % if any(ismember(fu.indeterminates.str,u.str))
      
        B = casos.PD(subs(fu,[x;u],[x0;u0]));

    % else
    %     B = full(fu);
    % end
else
    B = [];
end


end